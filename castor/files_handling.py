#! /usr/bin/env python3

import os
import re
import warnings

from astropy.io import fits
from tqdm import tqdm
import dateutil.parser
import numpy as np

FLOAT_DTYPE = np.float32

def get_package_data(path):
    package_root = os.path.abspath(os.path.dirname(__file__))
    return os.path.join(package_root, 'data', path)

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
