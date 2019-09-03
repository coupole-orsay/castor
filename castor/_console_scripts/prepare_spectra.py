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
        help=('Rotation angle of the spectrum, in degrees. '
              'Determined automatically if not specified.'))
    parser.add_argument(
        '--grating-period', type=int,
        help=('Period of the grating, in lines/mm. '
              'If not specified, --calib-points must be specified.'
            ))
    parser.add_argument(
        '--calib-points', type=str,
        help=('2-columns text file containing indices and wavelength of '
              'at least two wavelength calibration points. '
              'If not specified, --grating-period must be passed.'
            ))
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

    got_grating_period = (args.grating_period is not None)
    got_calib_points = (args.calib_points is not None)
    if not (got_grating_period or got_calib_points):
        parser.error('either --grating-period or --calib-points  are required')

    # output filenames
    args.sci_cube_r = os.path.join(args.target_name, 'cube_prepared_r.fits')
    args.master_calib = os.path.join(args.target_name, 'master_calib.fits')
    args.master_calib_r = os.path.join(args.target_name, 'master_calib_r.fits')

    return args

def main():
    args = get_parsed_args()

    if not args.overwrite:
        output_files = (
            args.sci_cube_r,
            args.master_calib,
            args.master_calib_r,
            )
        for output_file in output_files:
            if os.path.exists(output_file):
                msg = "output file '{}' exists, use -O to overwrite outputs"
                raise OSError(msg.format(output_file))

    master_calib = preparation.create_master(args.calib_path)
    # TODO: reduce size of master_calib to speed up RT
    # maybe this:
    import papy.num
    master_calib_small = papy.num.rebin(master_calib, (4, 4), cut_to_bin=True)
    # end TODO
    angle = spectroscopy.find_spectrum_orientation(master_calib_small)
    master_calib_r = sndi.rotate(master_calib, - angle)
    files_handling.save_fits(
        master_calib, args.master_calib, overwrite=args.overwrite)
    files_handling.save_fits(
        master_calib_r, args.master_calib_r, overwrite=args.overwrite)

    sci_cube_hdulist = fits.open(args.sci_cube)
    sci_cube = sci_cube_hdulist[args.hdu].data

    sci_cube_r = sndi.rotate(sci_cube, - angle, axes=(1, 2))

    sci_cube_hdulist[args.hdu].data = sci_cube_r
    sci_cube_hdulist.writeto(args.sci_cube_r, overwrite=args.overwrite)


if __name__ == '__main__':
    main()
