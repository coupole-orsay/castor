#!/usr/bin/env python3

import argparse
import os

from castor import files_handling, preparation

def get_parsed_args():
    parser = argparse.ArgumentParser(
        description='Prepare a series of images.',
        )
    parser.add_argument(
        'target_name', type=str,
        help='Name of the target')
    parser.add_argument(
        '--sci-path', type=str,
        help='Directory containing the science FITS. Default: {target_name}/sci')
    parser.add_argument(
        '--sci-dark-path', type=str,
        help='Directory containing the science dark FITS. Default: {target_name}/sci_dark')
    parser.add_argument(
        '--flat-path', type=str,
        help='Directory containing the flat FITS. Default: {target_name}/flat')
    parser.add_argument(
        '--flat-dark-path', type=str,
        help='Directory containing the flat dark FITS. Default: {target_name}/flat_dark')
    parser.add_argument(
        '-o', '--output', type=str,
        help='Output cube FITS. Default: {target_name}/cube_prepared.fits')
    parser.add_argument(
        '-O', '--overwrite', action='store_true',
        help='Overwrite output if it already exists.')

    args = parser.parse_args()

    # set defaults
    if args.sci_path is None:
        args.sci_path = os.path.join(args.target_name, 'sci')
    if args.flat_path is None:
        args.flat_path = os.path.join(args.target_name, 'flat')
    if args.sci_dark_path is None:
        args.sci_dark_path = os.path.join(args.target_name, 'sci_dark')
    if args.flat_dark_path is None:
        args.flat_dark_path = os.path.join(args.target_name, 'flat_dark')
    if args.output is None:
        args.output = os.path.join(args.target_name, 'cube_prepared.fits')

    return args

def main():
    args = get_parsed_args()

    if os.path.exists(args.output) and not args.overwrite:
        msg = "output file '{}' exists, use -O to overwrite it"
        raise OSError(msg.format(args.output))

    files_handling.compute_and_save(
        args.output, preparation.prepare,
        args.sci_path, args.sci_dark_path,
        args.flat_path, args.flat_dark_path,
        overwrite=args.overwrite,
        )

if __name__ == '__main__':
    main()
