''' Codes pour l’ASTronomie à ORsay

Collections of tools for the astronomy cupola in Orsay.

Submodules
==========
alignment :
    Tools to align (correct translation) and register (align + correct rotation
    and scaling) series of images.
photometry :
    Tools to extract sources from images using SEP.
spectroscopy :
    Tools to handle spectra acquired with the LHIRES III spectrometer.
files_handling :
    Helper functions to handle files, containing useful functions to load
    series of FITS.

In addition, the submodule mount_alignment (imported with `import
castor.mount_alignment`) contains tools to help align an equatorial mount.

Command line interfaces
=======================
These tools are mainly intended to be used through command-line scripts, which
are installed with pip (see README.md) :

castor_prepare :
    Prepare a series of images by correcting the dark and the flat.
castor_align :
    Align a series of images (not suitable for spectral data).
castor_rotate_spectra :
    Rotate spectral images so that the spectral and spatial (along the slit)
    dimensions of the 2D spectra are aligned with the axes of the images.
castor_wavelength_calibration :
    Perform the wavelength calibration for rotated spectra.
castor_exoplanet_analysis :
    Full pipeline for plotting exoplanet transit curves from raw data. Performs
    data preparation, alignment, source selection (using a GUI for the user to
    select the exoplanet host and reference stars), and photometry.
castor_pointing_analysis :
    Analyse the pointing stability of an equatorial mount from a series of
    images. Plot mount periodic error and drift, as well field rotation.
'''

from . import alignment, photometry, preparation, spectroscopy
