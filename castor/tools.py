#! /usr/bin/env python3

import argparse
import os
import re
import warnings

from astropy.io import fits
from tqdm import tqdm
import astroalign
import sep
import cv2
import dateutil.parser
import numpy as np

FLOAT_DTYPE = np.float32
PACKAGE_ROOT = os.path.abspath(os.path.dirname(__file__))

# Resource handling

def get_resource(path):
    return os.path.join(PACKAGE_ROOT, 'data', path)

# File handling

def list_fits(directory):
    if directory is None:
        return []
    try:
        all_files = os.listdir(directory)
    except FileNotFoundError:
        return []
    reg = re.compile('(?i).+\.fits?$')
    l = [os.path.join(directory, f)
        for f in all_files
        if reg.match(f)]
    l = sorted(l)
    return l

def load_fits_headers(filenames, hdu=0):
    headers = []
    for i, filename in enumerate(tqdm(filenames, desc='Loading headers')):
        f = fits.open(filename)
        f = f[hdu]
        headers.append(f.header)
    return headers

def load_fits_data(path, hdu=0, timestamps_hdu=None,
        norm_to_exptime=True, norm_dtype=FLOAT_DTYPE):
    f = fits.open(path)
    data = f[hdu].data
    if norm_to_exptime:
        data = data.astype(norm_dtype)
        data /= f[hdu].header['EXPTIME']
    if timestamps_hdu is not None:
        timestamps = f[timestamps_hdu].data['DATE-OBS']
        timestamps = np.array([dateutil.parser.parse(ts) for ts in timestamps])
        return data, timestamps
    else:
        return data

def get_timestamps(filenames, hdu=0):
    headers = load_fits_headers(filenames, hdu=hdu)
    timestamps = [dateutil.parser.parse(h['DATE-OBS']) for h in headers]
    return np.array(timestamps)

def open_or_compute(filename, function, *args, save=True, **kwargs):
    ''' If filename exists, open it; if it doesn't, compute it using
    function(*args, **kwargs) and save it to filename. '''
    try:
        data, timestamps = load_fits_data(filename,
            norm_to_exptime=False, timestamps_hdu=1)
    except FileNotFoundError:
        data, timestamps = function(*args, **kwargs)
        if save:
            try:
                save_fits(data, filename, timestamps=timestamps)
            except Exception as e:
                os.remove(filename)
                msg = '{} occured while saving {}: {}'
                warnings.warn(msg.format(e.__class__.__name__, filename, e))
    return data, timestamps

def save_fits(data, filename, overwrite=False, timestamps=None):
    hdulist = fits.HDUList([fits.PrimaryHDU(data)])
    if timestamps is not None:
        col = fits.Column(
            name='DATE-OBS',
            format='29A',
            array=timestamps,
            )
        hdulist.append(fits.BinTableHDU.from_columns([col]))
    hdulist.writeto(filename, overwrite=overwrite)

# Data preparation

def create_master(filenames, default=None, dtype=FLOAT_DTYPE):
    ''' Create a master dark or a master flat.

    Parameters
    ==========
    filenames : str or list
        If list, it is assumed to be a list of FITS filenames.
        If str, it is assumed to be a directory containing several FITS files.
    default : (default: None)
        Value returned if filenames is an empty directory.

    Returns
    =======
    master : 2D ndarray
        The average of all FITS in filenames, normalised to their respective
        exposure time.
    '''
    if type(filenames) is str:
        fits_path = filenames
        filenames = list_fits(fits_path)
    if not filenames:
        if default is not None:
            msg = '{} is empty, using default value {}'
            msg = msg.format(filenames, default)
            warnings.warn(msg)
            return default
        else:
            raise ValueError('{} is empty'.format(filenames))
    master = np.zeros_like(load_fits_data(filenames[0]))
    for filename in filenames:
        master += load_fits_data(filename)
    master = master / len(filenames)
    return master

def reduction(sci_dir, sci_dark_dir, flat_dir, flat_dark_dir):
    master_sci_dark = create_master(sci_dark_dir, default=0)
    master_flat_dark = create_master(flat_dark_dir, default=0)
    master_flat = create_master(flat_dir, default=1)
    master_flat = master_flat - master_flat_dark
    master_flat /= np.mean(master_flat)

    sci_filenames = list_fits(sci_dir)
    if not sci_filenames:
        raise ValueError('{} is empty'.format(sci_dir))

    # sort images with DATE-OBS
    timestamps = get_timestamps(sci_filenames)
    sci_order = np.argsort(timestamps)
    timestamps = timestamps[sci_order]
    sci_filenames = np.array(sci_filenames)[sci_order]

    # open sci data
    n_files = len(sci_filenames)
    sample_data = load_fits_data(sci_filenames[0])
    sci_images = np.zeros((n_files, *sample_data.shape), dtype=sample_data.dtype)
    for i, sci_filename in enumerate(tqdm(sci_filenames, desc='Opening FITS')):
        sci = load_fits_data(sci_filename)
        sci_images[i] = (sci - master_sci_dark) / master_flat
    return sci_images, timestamps

# Alignment

def affine_transform(img, mat):
    img_min = np.nanmin(img)
    img_max = np.nanmax(img)
    img = (img - img_min) / (img_max - img_min)
    img_transformed = cv2.warpAffine(
        img, mat, img.T.shape,
        borderValue=np.nan,
        )
    img_transformed = img_transformed * (img_max - img_min) + img_min
    return img_transformed

def register_stars(images, timestamps, ref_img=None):
    ''' Register a field of stars in translation, rotation, and scaling.

    Parameters
    ==========
    images : ndarray of shape (N, ny, nx)
        cube containing the images to align.
    ref_img : ndarray of shape (ny, nx) or None (default: None)
        The reference image relatively to which all images should be
        aligned.
        If None, use the first input image.
    Returns
    =======
    registered_images: ndarray of shape (N, ny, nx)
        version of the input, with all images aligned with ref_img.
    '''

    # registered_images = np.full_like(images, np.nan)

    if ref_img is None:
        ref_img = images[0]
        first_im_is_ref = True
    else:
        first_im_is_ref = False

    ref_sources = astroalign._find_sources(ref_img)
    iterable = tqdm(images, desc='Aligning images', total=len(images))
    for i, img in enumerate(iterable):
        if i == 0 and first_im_is_ref:
            continue
        try:
            p, _ = astroalign.find_transform(img, ref_sources)
            mat = p.params[:-1]
        except Exception as e:
            warnings.warn('Image {}: {}'.format(i, e))
            mat = np.array([[1, 0, 0], [0, 1, 0]], dtype=float)
        images[i] = affine_transform(img, mat)

    return images, timestamps

# Photometry

def sep_extract(data, threshold=3):
    ''' Extract sources from an image using SEP.

    Parameters
    ==========
    data : 2d ndarray
        Image containing the sources
    threshold : float
        The threshold value for detection, in number of sigma.

    Returns
    =======
    sources : np.recarray
        A list of sources, as returned by sep.extract, and ordered by flux.
        See documentation of sep.extract for a description of the fields.
    '''
    if isinstance(data, np.ma.MaskedArray):
        image = data.filled(fill_value=np.median(data)).astype(np.float32)
    else:
        image = data.astype(np.float32)
    bkg = sep.Background(image)
    thresh = threshold * bkg.globalrms
    sources = sep.extract(image - bkg.back(), thresh)
    sources.sort(order='flux')
    # sources = sources.view(np.recarray)
    return sources

def sep_sources_coordinates(*args, **kwargs):
    ''' Extract sources coordinates from an image using SEP.

    Parameters
    ==========
    *args and **kwarg : passed to sep_extract
    '''
    sources = sep_extract(*args, **kwargs)
    coordinates = [list(xy) for xy in sources[['x', 'y']]]
    return coordinates

def find_closest_sources(catalog, coordinates):
    ''' Find the sources of a catalog closest to a set of coordinates

    Parameters
    ==========
    catalog : (m, ) ndarray
        A sources catalog returned by sep.extract()
    coordinates : (n, 2) ndarray
        A list of x and y coordinates for the sources to find.

    Returns
    =======
    filtered_catalog : (n, ) ndarray
        A subset of the input catalog containing the sources which are the
        closest to the input coordinates
    distances : (n, ) ndarray
        The distances of each returned source to the input coordinate.
    '''

    # compute distance between all input and catalog sources
    cat_coordinates = np.stack((catalog['x'], catalog['y'])).T
    cat_coordinates  = cat_coordinates.reshape(-1, 2, 1) # (m, 2, 1)
    coordinates = np.array(coordinates)
    coordinates = coordinates.T.reshape(1, 2, -1) # (1, 2, n)
    dist = np.sum((cat_coordinates - coordinates)**2, axis=1) # (m, n)
    dist = np.sqrt(dist)

    # find closest sources of catalog
    i_cat, i_input = np.where(dist == dist.min(axis=0))

    # retrieve coordinates and distance in the correct order
    sort = np.argsort(i_input)
    i_cat = i_cat[sort]
    i_input = i_input[sort]
    filtered_catalog = catalog[i_cat]
    dist = dist[(i_cat, i_input)]

    return filtered_catalog, dist
