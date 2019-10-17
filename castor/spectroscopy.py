#!/usr/bin/env python3

import numpy as np
import scipy.interpolate as sinterp
from scipy.signal import fftconvolve
import skimage.transform as skt
from scipy.optimize import curve_fit
from tqdm import tqdm
from copy import deepcopy

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
    '''Generate an array of the pixel index - wavelength correspondence,
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

def cross_correlation(img1, img2):
    '''Determine the shift along the x and y axis between the two input
    2D images, using a 2D FFT convolution.

    Parameters
    ==========
    img1 : 2D ndarray
        The reference 2D spectrum.
    img2 : 2D ndarray
        The 2D spectrum to align

    Returns
    =======
    x_dith : int
        The shift along the x-axis in pixels.
    y_dith : int
        The shift along the y-axis in pixels.
    '''
    ysize, xsize = img1.shape
    corr = fftconvolve(img1, img2[::-1, ::-1], mode='same')
    y_dith, x_dith = np.where(corr == np.max(corr))
    x_dith = int(x_dith) - xsize//2
    y_dith = int(y_dith) - ysize//2
    return x_dith, y_dith

def align_images(Img, img_ref=0, xaxis=True, yaxis=True, 
                 gauss=True, sigma=0.5):
    '''Align a series of image relatively to the img_ref,
    along the x and y-axis.

    Parmeters
    =========
    Img : 3D ndarray
        The array of 2D spectra (images) to align.
    img_ref : int (default: 0)
        The index of the image taken as a reference.
    xaxis : bool (default: True)
        Enable x-axis alignment.
    yaxis : bool (default: True)
        Enable y-axis alignment.
    gauss : bool (default: True)
        If True, the images are convolve with a 2D gaussian
        to reduce the effects of sharpy edges.
    sigma : float (default : 0.5)
        If gauss=True, the standard deviation of the use gaussian.

    Returns
    =======
    Img_align : 3D ndarray
        The aligned array of 2D spectra.
    '''
    N, Ysize, Xsize = Img.shape
    Img2 = deepcopy(Img)
    # Gaussian convolution
    if gauss:
        X, Y = np.linspace(-1, 1, Xsize), np.linspace(-1, 1, Ysize)
        G = gauss2d(X, Y, 1, 0, 0, sigma, sigma)
        for i in range(N):
            Img2 *= G
    # Shift determination using FFT convolution
    ref = Img2[img_ref]
    x_dith, y_dith = np.zeros(N, dtype=int), np.zeros(N, dtype=int)
    for i in tqdm(range(N)):
        x_dith[i], y_dith[i] = cross_correlation(ref, Img2[i])
    if not xaxis:
        x_dith[:] = 0
    if not yaxis:
        y_dith[:] = 0
    # Creation of a new set of images, with new dimensions
    x0 = np.max(np.abs(x_dith))
    y0 = np.max(np.abs(y_dith))
    Xnew = Xsize + 2*x0
    Ynew = Ysize + 2*y0
    Img_align = np.zeros((N, Ynew, Xnew))
    # Alignement
    for i in range(N):
        Img_align[i, y0+y_dith[i]:y0+Ysize+y_dith[i], x0+x_dith[i]:x0+Xsize+x_dith[i]] = Img[i]
    # Output
    return Img_align

def gauss2d(x, y, a, x0, y0, sigmax, sigmay):
    '''Compute a 2D gaussian on a grid generated from the
    input x and y arrays.

    G = a * exp( -(x-x0)**2/(2*sigmax**2) - (y-y0)**2/(2*sigmay**2) )

    Parameters
    ==========
    x : 1D ndarray
        The x-axis values.
    y : 1D ndarray
        The y-axis values.
    a : float
        The amplitude of the gaussian.
    x0 : float
        The x-coordinate of the gaussian center.
    y0 : float
        The y-coordinate of the gaussian center.
    sigmax : float
        The standard deviation along the x-axis.
    sigmay : float
        The standard deviation along the y-axis.

    Returns
    =======
    G : 2D ndarray
        The 2D gaussian computed of the (x, y) grid.
    '''
    X, Y = np.meshgrid(x, y)
    G = a * np.exp( -(X-x0)**2/(2*sigmax**2) - (Y-y0)**2/(2*sigmay**2) )
    return G
