#!/usr/bin/env python3

import warnings

from tqdm import tqdm
import numpy as np

from . import files_handling

def create_master(filenames, default=None):
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
        filenames = files_handling.list_fits(fits_path)
        empty_name = fits_path
    else:
        empty_name = filenames
    if not filenames:
        if default is not None:
            msg = '{} is empty, using default value {}'
            msg = msg.format(empty_name, default)
            warnings.warn(msg)
            return default
        else:
            raise ValueError('{} is empty'.format(empty_name))
    master = np.zeros_like(files_handling.load_fits_data(filenames[0]))
    for filename in filenames:
        master += files_handling.load_fits_data(filename)
    master = master / len(filenames)
    return master

def prepare(sci_dir, sci_dark_dir, flat_dir, flat_dark_dir):
    master_sci_dark = create_master(sci_dark_dir, default=0)
    master_flat_dark = create_master(flat_dark_dir, default=0)
    master_flat = create_master(flat_dir, default=1)
    master_flat = master_flat - master_flat_dark
    master_flat /= np.mean(master_flat)

    sci_filenames = files_handling.list_fits(sci_dir)
    if not sci_filenames:
        raise ValueError('{} is empty'.format(sci_dir))

    # sort images with DATE-OBS
    timestamps = files_handling.get_timestamps(sci_filenames)
    sci_order = np.argsort(timestamps)
    timestamps = timestamps[sci_order]
    sci_filenames = np.array(sci_filenames)[sci_order]

    # open sci data
    n_files = len(sci_filenames)
    sample_data = files_handling.load_fits_data(sci_filenames[0])
    sci_images = np.zeros((n_files, *sample_data.shape), dtype=sample_data.dtype)
    for i, sci_filename in enumerate(tqdm(sci_filenames, desc='Opening FITS')):
        sci = files_handling.load_fits_data(sci_filename)
        sci_images[i] = (sci - master_sci_dark) / master_flat
    return sci_images, timestamps
