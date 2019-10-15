#!/usr/bin/env python3

import os
import re
import warnings

from astropy.io import fits
from tqdm import tqdm
import dateutil.parser
import numpy as np

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
        norm_to_exptime=True, norm_dtype=np.float32):
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

def pass_timestamps(func):
    ''' Decorator to make a function transparently pass timestamps

    Used to make functions that don't use timestamps (eg.
    alignment.register_stars) work with open_or_compute.

    transforms:
        ret = func(a, b, ..., c=foo, d=bar)
    into:
        ret, timestamps = func(a, timestamps, b, ..., c=foo, d=bar)
    '''
    def func_that_passes_timestamps(a, timestamps, *args, **kwargs):
        return func(a, *args, **kwargs), timestamps
    return func_that_passes_timestamps

def compute_and_save(filename, function, *args, overwrite=False, **kwargs):
    data, timestamps = function(*args, **kwargs)
    save_fits(data, filename, timestamps=timestamps, overwrite=overwrite)

def open_or_compute(filename, function, *args,
                    save=True, use_timestamps=True, **kwargs):
    ''' If filename exists, open it; if it doesn't, compute it using
    function(*args, **kwargs) and save it to filename.

    Parameters
    ==========
    filename : str
        FITS to open if it exists, or where to save the computed data
        if `save=True` is passed
    function : callable
        Callable which accepts *args and **kwargs, and returns either
        `data`, or `data, timestamps`
    save : bool (default: True)
        Whether to save the computed data.
    use_timestamps : bool (default: True)
        If True, assume that `function` returns timestamps.
    *args, **kwargs :
        Passed to `function`

    Returns
    =======
    data : ndarray
        The loaded or computed data
    timestamps : ndarray or None
        The timestamps computed or loaded from the FITS.
        None if `use_timestamps=False`.
    '''
    timestamps = None
    if use_timestamps:
        timestamps_hdu = 1
    else:
        timestamps_hdu = None
    if os.path.exists(filename):
        data = load_fits_data(
            filename,
            norm_to_exptime=False,
            timestamps_hdu=timestamps_hdu,
            )
        if use_timestamps:
            data, timestamps = data
    else:
        data = function(*args, **kwargs)
        if use_timestamps:
            data, timestamps = data
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
