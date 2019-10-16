#!/usr/bin/env python3

import argparse
import os
import numpy as np
from astropy.io import fits
from castor import spectroscopy

def get_parsed_args():
    parser = argparse.ArgumentParser(
        description='Align a cube of 2D spectra.',
        )
    parser.add_argument(
        'target_name', type=str,
        help='Name of the target')
    parser.add_argument(
        '--sci-cube-rotated', type=str,
        help=('Science data prepared with castor_prepare, '
              'and rotated with castor_rotate_spectra.'
              ' Default: {target_name}/cube_prepared_rotated.fits'))
    parser.add_argument(
        '--output', type=str,
        help='Directory where output files are saved. Default: {target_name}/cube_prepared_rotated_aligned.fits')
    parser.add_argument(
        '-O', '--overwrite', action='store_true',
        help='Overwrite outputs if they already exist.')
    parser.add_argument(
        '--hdu', type=int, default=0,
        help='HDU containing the cube to align. Default: 0')
    args = parser.parse_args()

    # default inputs and outputs
    if not args.sci_cube_rotated:
        args.sci_cube_rotated = os.path.join(args.target_name, 'cube_prepared_rotated.fits')
    if not args.output:
        args.output = os.path.join(args.target_name, 'cube_prepared_rotated_aligned.fits')

    return args

def main():
    args = get_parsed_args()

    if not args.overwrite:
        if os.path.exists(args.output):
            msg = "output file '{}' exists, use -O to overwrite it"
            raise OSError(msg.format(args.output))
    
    hdulist = fits.open(args.sci_cube_rotated)
    cube = hdulist[0].data
    cube_aligned = spectroscopy.align_images(cube, args.hdu, xaxis=False, yaxis=True)

    hdulist[args.hdu].data = cube_aligned
    hdulist.writeto(args.output, overwrite=args.overwrite)

if __name__ == '__main__':
    main()
