#!/usr/bin/env python3

import numpy as np
import scipy.interpolate as sinterp
import skimage.transform as skt
from scipy.optimize import curve_fit

def resize_for_radon(spectrum, min_size=150):
    ''' Resize 2D spectrum to speed up skt.radon

    Parameters
    ==========
    spectrum : 2D ndarray
        The 2D spectrum to resize.
    min_size : int (default: 150)
        Minimum size for the smallest axis of the resizedn array.
        (Note: quick testing indicates that using a spectrum scaled down to
        100x133 yields the same angle as the original 600x800 spectrum. Using a
        150x200 spectrum seems a good time/robustness tradeof, as it takes
        approximately 5s to compute skt.radon with a 600x800 array, and 1s with
        a 150x200 array.)

    Returns
    =======
    new_spectrum : 2D ndarray
        A sized-down version of the input spectrum.
    '''
    shape = np.array(spectrum.shape)
    smallest_axis_size = np.min(shape)
    if smallest_axis_size <= min_size:
        return spectrum
    # compute new shape
    new_shape = shape / smallest_axis_size * min_size
    new_shape = np.round(new_shape).astype(int)
    # build x and y coordinates for the input spectrum
    ny, nx = shape
    new_ny, new_nx = new_shape
    y = np.arange(ny)
    x = np.arange(nx)
    # buidl new_x and new_y coordinates for the resized spectrum
    new_y = np.linspace(y[0], y[-1], new_ny)
    new_x = np.linspace(x[0], x[-1], new_nx)
    # interpolate the new spectrum
    interp_spectrum = sinterp.interp2d(x, y, spectrum)
    new_spectrum = interp_spectrum(new_x, new_y)
    return new_spectrum

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
    spectrum_small = resize_for_radon(spectrum)
    angles = np.arange(0, 180, angle_step)
    # Radon transform: axis 0: displacement; axis1: angle
    spectrum_rt = skt.radon(spectrum_small, angles, circle=False)
    # maximum of the RT across all displacements:
    spectrum_rt_max = np.max(spectrum_rt, axis=0)
    # for a spectrum compose dof straight emission lines, the global
    # maximum of the RT gives the orientation of the lines.
    angle = angles[np.argmax(spectrum_rt_max)]
    return angle

def calib_wavelength_array(calib_pts, Nlam):
    '''Generate a array of the pixel index - wavelength correspondence,
    from an linear fit of some (pixel_index, associated wavelenth) tuple.

    Parameters
    ==========
    calib_pts : 2D ndarray
        A 2D array containing the pixel index and the associated
        wavelength (at least 2 calibration points are required).
    Nlam : int
        The total number of pixels along the wavelength axis.

    Returns
    =======
    calib_array : 2D ndarray
        A 2D array containing for each pixel along the wavelength axis
        the assoicated wavelength.
    '''
    # Initialization
    px_array = np.arange(Nlam)
    # Linear fitting
    f_lin = lambda x, a, b : a*x + b
    a, b = curve_fit(f_lin, calib_pts[:,0], calib_pts[:,1])[0]
    lam_array = a * px_array + b
    # Output
    calib_array = np.array([px_array, lam_array]).T
    return calib_array
