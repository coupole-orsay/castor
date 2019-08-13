#!/usr/bin/env python3

import warnings

from tqdm import tqdm
import astroalign
import cv2
import numpy as np

def affine_transform(img, mat):
    img_min = np.nanmin(img)
    img_max = np.nanmax(img)
    img = (img - img_min) / (img_max - img_min)
    img_transformed = cv2.warpAffine(
        img, mat, img.T.shape,
        borderValue=np.nan,
        )
    img_transformed = img_transformed * (img_max - img_min) + img_min
    return img_transformed

def register_stars(images, ref_img=None):
    ''' Register a field of stars in translation, rotation, and scaling.

    Parameters
    ==========
    images : ndarray of shape (N, ny, nx)
        cube containing the images to align.
    ref_img : ndarray of shape (ny, nx) or None (default: None)
        The reference image relatively to which all images should be
        aligned.
        If None, use the first input image.
    Returns
    =======
    registered_images: ndarray of shape (N, ny, nx)
        version of the input, with all images aligned with ref_img.
    '''

    # registered_images = np.full_like(images, np.nan)

    if ref_img is None:
        ref_img = images[0]
        first_im_is_ref = True
    else:
        first_im_is_ref = False

    ref_sources = astroalign._find_sources(ref_img)
    iterable = tqdm(images, desc='Aligning images', total=len(images))
    for i, img in enumerate(iterable):
        if i == 0 and first_im_is_ref:
            continue
        try:
            p, _ = astroalign.find_transform(img, ref_sources)
            mat = p.params[:-1]
        except Exception as e:
            warnings.warn('Image {}: {}'.format(i, e))
            mat = np.array([[1, 0, 0], [0, 1, 0]], dtype=float)
        images[i] = affine_transform(img, mat)

    return images
