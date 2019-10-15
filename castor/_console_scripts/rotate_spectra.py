#!/usr/bin/env python3

import argparse
import os

from astropy.io import fits
import scipy.ndimage as sndi

from castor import files_handling, preparation, spectroscopy

def get_parsed_args():
    parser = argparse.ArgumentParser(
        description='Prepare a series of images.',
        )
    parser.add_argument(
        'target_name', type=str,
        help='Name of the target')
    parser.add_argument(
        '--spectrum-rotation', type=float,
        default=None,
        help=('Rotation angle of the spectrum, in degrees. '
              'Determined automatically if not specified.'))
    parser.add_argument(
        '--sci-cube', type=str,
        help=('Science data prepared with castor_prepare.'
              ' Default: {target_name}/cube_prepared.fits'))
    parser.add_argument(
        '--calib-path', type=str,
        help='Directory containing the spectral calibration FITS. Default: {target_name}/calib')
    parser.add_argument(
        '--output-path', type=str,
        help='Directory where output files are saved. Default: {name}/')
    parser.add_argument(
        '-O', '--overwrite', action='store_true',
        help='Overwrite outputs if they already exist.')
    parser.add_argument(
        '--hdu', type=int, default=0,
        help='HDU containing the cube to align. Default: 0')

    args = parser.parse_args()

    # set defaults
    if not args.calib_path:
        args.calib_path = os.path.join(args.target_name, 'calib')
    if not args.sci_cube:
        args.sci_cube = os.path.join(args.target_name, 'cube_prepared.fits')
    if not args.output_path:
        args.output_path = args.target_name

    # output filenames
    args.master_calib = os.path.join(args.target_name, 'master_calib.fits')
    args.master_calib_rotated = os.path.join(args.target_name, 'master_calib_rotated.fits')
    args.sci_cube_rotated = os.path.join(args.target_name, 'cube_prepared_rotated.fits')

    return args

def main():
    args = get_parsed_args()

    if not args.overwrite:
        output_files = (
            args.sci_cube_rotated,
            args.master_calib_rotated,
            )
        for output_file in output_files:
            if os.path.exists(output_file):
                msg = "output file '{}' exists, use -O to overwrite outputs"
                raise OSError(msg.format(output_file))

    # make master calib -------------------------------------------------------
    master_calib, _ = files_handling.open_or_compute(
        args.master_calib,
        preparation.create_master,
        args.calib_path,
        use_timestamps=False,
        )

    # set or find spectrum rotation angle -------------------------------------
    if args.spectrum_rotation is not None:
        angle = args.spectrum_rotation
        print('Using provided spectrum rotation: {:.2f}°'.format(angle))
    else:
        angle = spectroscopy.find_spectrum_orientation(master_calib)
        print('Using computed spectrum rotation: {:.2f}°'.format(angle))

    # rotate calib and science spectra ----------------------------------------

    master_calib_rotated = sndi.rotate(master_calib, - angle)
    files_handling.save_fits(
        master_calib_rotated, args.master_calib_rotated, overwrite=args.overwrite)

    sci_cube_hdulist = fits.open(args.sci_cube)
    sci_cube = sci_cube_hdulist[args.hdu].data

    sci_cube_rotated = sndi.rotate(sci_cube, - angle, axes=(1, 2))

    sci_cube_hdulist[args.hdu].data = sci_cube_rotated
    sci_cube_hdulist.writeto(args.sci_cube_rotated, overwrite=args.overwrite)


if __name__ == '__main__':
    main()
