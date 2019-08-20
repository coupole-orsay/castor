#!/usr/bin/env python3

import argparse
import os
import warnings

from astropy import units as u
from astropy.io import fits
from tqdm import tqdm
import astroalign
import numpy as np
import skimage.transform as skt
import yaml

from castor import files_handling, preparation, alignment, photometry

def get_parsed_args():
    ''' Script argument parser '''
    parser = argparse.ArgumentParser(
        description='Telescope alignment analysis.')
    parser.add_argument(
        'name', type=str,
        help='Name of the alignment measurement')
    parser.add_argument(
        '--track-path', type=str,
        help='Directory containing the tracking FITS. Default: {name}/track')
    parser.add_argument(
        '--sample-track-fits', type=str,
        help='Sample track FITS within the track path. Default: first FITS in the track path.')
    parser.add_argument(
        '--track-dark-path', type=str,
        help='Directory containing the tracking dark FITS. Default: {name}/track_dark')
    parser.add_argument(
        '--flat-path', type=str,
        help='Directory containing the flat FITS. Default: {name}/flat')
    parser.add_argument(
        '--flat-dark-path', type=str,
        help='Directory containing the flat dark FITS. Default: {name}/flat_dark')
    parser.add_argument(
        '--output-path', type=str,
        help='Directory where output files are saved. Default: {name}/')
    parser.add_argument(
        '--sep-threshold', type=float, default=20,
        help=('Source extraction threshold, passed to sep.extract(). '
              'https://sep.readthedocs.io/en/latest/api/sep.extract.html '
              'Default: 20'))

    args = parser.parse_args()

    if args.track_path is None:
        args.track_path = os.path.join(args.name, 'track')
    if args.sample_track_fits is None:
        track_files = files_handling.list_fits(args.track_path)
        if not track_files:
            raise ValueError('No track files found in ' + args.track_path)
        args.sample_track_fits = track_files[0]
    if args.flat_path is None:
        args.flat_path = os.path.join(args.name, 'flat')
    if args.track_dark_path is None:
        args.track_dark_path = os.path.join(args.name, 'track_dark')
    if args.flat_dark_path is None:
        args.flat_dark_path = os.path.join(args.name, 'flat_dark')
    if args.output_path is None:
        args.output_path = args.name

    return args

def find_transform(sources_1, sources_2):
    try:
        p, _ = astroalign.find_affine_transform(sources_1, sources_2)
    except Exception as e:
        warnings.warn(str(e))
        p = skt.SimilarityTransform(
            scale=np.nan, rotation=np.nan, translation=np.nan)
    return p

def main():
    args = get_parsed_args()

    # Instruments properties
    with open(files_handling.get_package_data('instruments.yml')) as f:
        instruments = yaml.load(f)
    focale = u.Quantity(instruments['telescopes']['c14']['focale length'])
    px_size = instruments['cameras']['atik']['pixel size']
    px_size = u.Quantity([u.Quantity(s.replace('µm', 'um')) for s in px_size])

    # Data properties
    sample_fits = os.path.join(args.sample_track_fits)
    hdu = fits.open(sample_fits)[0]
    binning = np.array((hdu.header['XBINNING'], hdu.header['YBINNING']))
    binned_px_size = px_size * binning
    px_angular_size = np.arctan(binned_px_size / focale).to('arcsec')

    # Open data
    cube_path = os.path.join(args.output_path, 'cube_prepared.fits')
    images, timestamps = files_handling.open_or_compute(
        cube_path, preparation.prepare,
        args.track_path, args.track_dark_path,
        args.flat_path, args.flat_dark_path,
        # save=False,
        )

    # Align images ------------------------------------------------------------

    # Align using a reference image
    ref_img = images[0]
    ref_sources = photometry.sep_sources_coordinates(ref_img, threshold=args.sep_threshold)
    iterable = tqdm(images, desc='Aligning images', total=len(images))
    transforms = []
    for i, img in enumerate(iterable):
        try:
            sources = photometry.sep_sources_coordinates(img, threshold=args.sep_threshold)
            p, _ = astroalign.find_transform(sources, ref_sources)
        except Exception as e:
            warnings.warn('Image {}: {}'.format(i, e))
            p = skt.SimilarityTransform(
                scale=np.nan, rotation=np.nan, translation=np.nan)
        transforms.append(p)

    scale = np.array([p.scale for p in transforms])
    rotation = np.rad2deg(np.array([p.rotation for p in transforms]))
    translation = np.array([p.translation for p in transforms])
    translation = translation * px_angular_size.to('arcsec').value
    tx, ty = translation.T

    # plots -------------------------------------------------------------------
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    plt.ioff()

    style = dict(
        marker='.',
        linestyle='',
        )
    ref_date = timestamps[0]
    x_label = 'Time since {} [h]'.format(ref_date.isoformat(timespec='seconds'))
    rel_time = [(ts - ref_date).total_seconds() for ts in timestamps]
    rel_time = np.array(rel_time) / 3600

    plt.clf()
    plt.plot(rel_time, scale, **style)
    plt.xlabel(x_label)
    plt.ylabel('Scale [px]')
    plt.savefig(os.path.join(args.output_path, 'scale.pdf'))

    plt.clf()
    plt.plot(rel_time, rotation, **style)
    plt.xlabel(x_label)
    plt.ylabel('Rotation [°]')
    plt.savefig(os.path.join(args.output_path, 'rotation.pdf'))

    plt.clf()
    plt.plot(rel_time, tx, label='X', **style)
    plt.plot(rel_time, ty, label='Y', **style)
    plt.legend()
    plt.xlabel(x_label)
    plt.ylabel('Translation [arcsec]')
    plt.savefig(os.path.join(args.output_path, 'translation.pdf'))

    from scipy.optimize import curve_fit
    plt.clf()
    plt.legend()
    plt.xlabel(x_label)
    plt.ylabel('Translation [arcsec]')
    y_fill = np.linspace(*plt.ylim(), 100)
    fit_segments = [
        ('X', tx,
            lambda t, a, b, A, T, φ: a*t + b + A*np.sin(2*np.pi*t/T+φ),
            'X(t) = {:.2f} t + {:.2f} + {:.2f} sin(2 π t / {:.4f} + {:.2f})',
            (1, 0, 1, .083, 0)),
        ('Y', ty,
            lambda t, a, b: a*t + b,
            'Y(t) = {:.2f} t + {:.2f}',
            None),
        ]
    popts, pcovs = [], []
    fit_xdata = rel_time
    for varname, fx, func, expr, p0 in fit_segments:
        fit_ydata = fx
        m = ~np.isnan(fit_ydata)
        popt, pcov = curve_fit(func, fit_xdata[m], fit_ydata[m], p0=p0)
        popts.append(popt)
        pcovs.append(pcov)
        equation = expr.format(*popt)
        print(equation)
        plt.plot(
            rel_time, fit_ydata,
            label=equation,
            **style)
        plt.plot(
            fit_xdata[m], func(fit_xdata, *popt)[m],
            color='#808080',
            )
    plt.legend()
    plt.savefig(os.path.join(args.output_path, 'translation_fit.pdf'))

if __name__ == '__main__':
    main()
