#!/usr/bin/env python3

import sep
import numpy as np

def sep_extract(data, threshold=3):
    ''' Extract sources from an image using SEP.

    Parameters
    ==========
    data : 2d ndarray
        Image containing the sources
    threshold : float
        The threshold value for detection, in number of sigma.

    Returns
    =======
    sources : np.recarray
        A list of sources, as returned by sep.extract, and ordered by flux.
        See documentation of sep.extract for a description of the fields.
    '''
    if isinstance(data, np.ma.MaskedArray):
        image = data.filled(fill_value=np.median(data)).astype(np.float32)
    else:
        image = data.astype(np.float32)
    bkg = sep.Background(image)
    thresh = threshold * bkg.globalrms
    sources = sep.extract(image - bkg.back(), thresh)
    sources.sort(order='flux')
    # sources = sources.view(np.recarray)
    return sources

def sep_sources_coordinates(*args, **kwargs):
    ''' Extract sources coordinates from an image using SEP.

    Parameters
    ==========
    *args and **kwarg : passed to sep_extract
    '''
    sources = sep_extract(*args, **kwargs)
    coordinates = [list(xy) for xy in sources[['x', 'y']]]
    return coordinates

def find_closest_sources(catalog, coordinates):
    ''' Find the sources of a catalog closest to a set of coordinates

    Parameters
    ==========
    catalog : (m, ) ndarray
        A sources catalog returned by sep.extract()
    coordinates : (n, 2) ndarray
        A list of x and y coordinates for the sources to find.

    Returns
    =======
    filtered_catalog : (n, ) ndarray
        A subset of the input catalog containing the sources which are the
        closest to the input coordinates
    distances : (n, ) ndarray
        The distances of each returned source to the input coordinate.
    '''

    # compute distance between all input and catalog sources
    cat_coordinates = np.stack((catalog['x'], catalog['y'])).T
    cat_coordinates  = cat_coordinates.reshape(-1, 2, 1) # (m, 2, 1)
    coordinates = np.array(coordinates)
    coordinates = coordinates.T.reshape(1, 2, -1) # (1, 2, n)
    dist = np.sum((cat_coordinates - coordinates)**2, axis=1) # (m, n)
    dist = np.sqrt(dist)

    # find closest sources of catalog
    i_cat, i_input = np.where(dist == dist.min(axis=0))

    # retrieve coordinates and distance in the correct order
    sort = np.argsort(i_input)
    i_cat = i_cat[sort]
    i_input = i_input[sort]
    filtered_catalog = catalog[i_cat]
    dist = dist[(i_cat, i_input)]

    return filtered_catalog, dist
