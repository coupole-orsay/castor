#!/usr/bin/env python3

import numpy as np
import skimage.transform as skt

def find_spectrum_orientation(spectrum, angle_step=.25):
    ''' Determine the orientation of a 2D spectrum

    This uses a Radon transform to find the orientation of lines in an image.

    Parameters
    ==========
    spectrum : 2D ndarray
        A 2D spectrum, where dimensions are the position along the slit and the
        wavelength, and are not aligned with the array axes.
        The spectrum must contain emission lines that occupy the larger part of
        the slit height (eg. a spectrum of the ArNe calibration lamp).
    angle_step : float (default: .25)
        The angle step (in degrees) used when computing the Radon transform.
        This is roughly the accuracy of the resulting angle.

    Returns
    =======
    angle : float
        The rotation between the emission lines and the vertical axis of the
        input array.
        (Rotating `spectrum` by `- angle` would return an array where axis 0 is
        the position along the slit, and axis 1 the wavelength dimension.)
    '''
    angles = np.arange(0, 180, angle_step)
    # Radon transform: axis 0: displacement; axis1: angle
    spectrum_rt = skt.radon(spectrum, angles, circle=False)
    # maximum of the RT across all displacements:
    spectrum_rt_max = np.max(spectrum_rt, axis=0)
    # for a spectrum compose dof straight emission lines, the global
    # maximum of the RT gives the orientation of the lines.
    angle = angles[np.argmax(spectrum_rt_max)]
    return angle
