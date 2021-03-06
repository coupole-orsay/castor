#!/usr/bin/env python3

import argparse
import os
import numpy as np
from astropy.io import fits
from castor import spectroscopy

def get_parsed_args():
    parser = argparse.ArgumentParser(
        description='Compute the wavelength calibration.',
        )
    parser.add_argument(
        'target_name', type=str,
        help='Name of the target')
    parser.add_argument(
        '--calib-points', type=str,
        help=('2-columns text file containing indices and wavelength of '
              'at least two wavelength calibration points. '
              'Each line contains the indice followed by a space '
              'and the corresponding wavelength. '
              'Default: {name}/wavelength_calib_points.txt'))
    parser.add_argument(
        '--start-index', type=int, default=0,
        help=('1st index of the wavelength array. '
              'Must be 0 or 1. Default: 0'))
    parser.add_argument(
        '--master-calib-rotated', type=str,
        help=('Rotated wavelength calibration image generated by rotate_spectra'
              ' Default: {target_name}/master_calib_rotated.fits'))
    parser.add_argument(
        '--output', type=str,
        help='Directory where output files are saved. Default: {name}/wavelength_array.txt')
    parser.add_argument(
        '-O', '--overwrite', action='store_true',
        help='Overwrite outputs if they already exist.')
    args = parser.parse_args()

    # default inputs and outputs
    if not args.calib_points:
        args.calib_points = os.path.join(args.target_name, 'wavelength_calib_points.txt')
    if not args.master_calib_rotated:
        args.master_calib_rotated = os.path.join(args.target_name, 'master_calib_rotated.fits')
    if not args.output:
        args.output = os.path.join(args.target_name, 'wavelength_array.txt')

    return args

def main():
    args = get_parsed_args()

    if not args.overwrite:
        if os.path.exists(args.output):
            msg = "output file '{}' exists, use -O to overwrite it"
            raise OSError(msg.format(args.output))
    
    f = fits.open(args.master_calib_rotated)
    data = f[0].data
    Nlam = data.shape[1]
    calib_pts = np.loadtxt(args.calib_points)
    calib_pts[:,0] -= args.start_index  # Re-index the array to start at 0
    calib_array = spectroscopy.calib_wavelength_array(calib_pts, Nlam)
    np.savetxt(args.output, calib_array)

if __name__ == '__main__':
    main()
