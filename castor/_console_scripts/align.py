#! /usr/bin/env python3

import argparse
import os

from castor import files_handling, alignment

def get_parsed_args():
    parser = argparse.ArgumentParser(
        description='Align a cube of images',
        )
    parser.add_argument(
        'target_name', type=str,
        nargs='?',
        help='Name of the target')
    parser.add_argument(
        '-i', '--input', type=str,
        help='Input cube FITS. Default: {target_name}/cube_prepared.fits')
    parser.add_argument(
        '-o', '--output', type=str,
        help='Output cube FITS. Default: {target_name}/cube_aligned.fits')
    parser.add_argument(
        '-O', '--overwrite', action='store_true',
        help='Overwrite output if it already exists.')

    args = parser.parse_args()

    # set defaults
    got_target_name = (args.target_name is not None)
    got_input_and_output = ((args.input is not None) and (args.output is not None))
    if not (got_target_name or got_input_and_output):
        parser.error('either target_name, or -i INPUT and -o OUTPUT are required')
    if args.input is None:
        args.input = os.path.join(args.target_name, 'cube_prepared.fits')
    if args.output is None:
        args.output = os.path.join(args.target_name, 'cube_aligned.fits')

    return args

def main():
    args = get_parsed_args()

    if os.path.exists(args.output) and not args.overwrite:
        msg = "output file '{}' exists, use -O to overwrite it"
        raise OSError(msg.format(args.output))

    images, timestamps = files_handling.load_fits_data(
        args.input, norm_to_exptime=False, timestamps_hdu=1)
    files_handling.compute_and_save(
        args.output,
        files_handling.pass_timestamps(alignment.register_stars),
        images, timestamps, overwrite=args.overwrite)

if __name__ == '__main__':
    main()
