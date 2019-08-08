#!/usr/bin/env python3

import numpy as np
from numpy import sin, cos
from astropy import units as u
from astropy import constants as c

Ωsid = u.Quantity(360, 'deg') / u.Quantity(24, 'h') * 366.2422 / 365.2422

def r2d(ang):
    try:
        if ang.to('').unit is u.dimensionless_unscaled:
            ang = ang.to('').value
    except (AttributeError, u.UnitconversionError):
        pass
    ang = u.Quantity(ang, 'rad').to('deg')
    return ang

def rotation_to_screw_turns(angle_variation, axis):
    if axis.lower() in ('az', 'azimuth', 'lon', 'longitude'):
        a = u.Quantity(3.400, 'cycle / deg')
    elif axis.lower() in ('alt', 'altitude', 'lat', 'latitude'):
        a = u.Quantity(3.707, 'cycle / deg')
    else:
        raise ValueError('Unknown axis: {}'.format(axis))

    return (a * angle_variation).to('cycle')

def drifts_to_adjustments(ωδ, ωH, δ, H, φ):
    ''' Get the alignment corrections to apply to a mount given measured
    drifts.

    Parameters
    ==========
    ωδ : declination drift
    ωH : right ascension drift
    δ : declination of the tracked star
    H : hour angle of the tracked star
    φ : latitude of observation location

    Returns
    =======
    Δφ : correction to the altitude axis
    ΔA : correction to the azimuth axis
    '''
    a = 1 / (Ωsid * sin(δ) * (sin(H)**2 + cos(H)**2 * cos(φ)))
    Δφ = a * (- sin(H) * sin(δ) * ωδ + cos(H) * cos(φ) * ωH)
    ΔA = a * (cos(H) * sin(δ) * ωδ + sin(H) * ωH)
    return r2d(Δφ), r2d(ΔA)

def bigourdan(ωδ, φ, where):
    if where.lower() in ('s', 'south'):
        ΔA = ωδ / (Ωsid * cos(φ))
        return r2d(ΔA)
    elif where.lower() in ('e', 'east'):
        Δφ = + ωδ / Ωsid
        return r2d(Δφ)
    elif where.lower() in ('w', 'west'):
        Δφ = - ωδ / Ωsid
        return r2d(Δφ)
    else:
        raise ValueError('unknown_location')

if __name__ == '__main__':

    # Coupole 470
    φ = u.Quantity(48.70272, 'deg')
    θ = u.Quantity(2.18270, 'deg')

    ωδ = u.Quantity(726.87, 'arcsec / h')
    ωH = u.Quantity(-31.78, 'arcsec / h')

    # # altair
    # δ = u.Quantity(9, 'deg')
    # H = u.Quantity(3.25, 'hourangle')

    # # deneb
    # δ = u.Quantity(45, 'deg')
    # H = u.Quantity(2.5, 'hourangle')

    # 59 aql 2018-11-30 19:20 CET
    δ = u.Quantity(8.5, 'deg')
    H = u.Quantity(3.25, 'hourangle')

    Δφ, ΔA = drifts_to_adjustments(ωδ, ωH, δ, H, φ)
    print('Alt  /  Az')
    print('{:.3g} {:.3g}'.format(Δφ, ΔA))
    print('{:.3g} {:.3g}'.format(
        rotation_to_screw_turns(Δφ, 'alt'), 
        rotation_to_screw_turns(ΔA, 'az'),
        ))

    # ωδ = u.Quantity(726.87, 'arcsec / h')
    # where = 'west'
    # Δ = bigourdan(ωδ, φ, where)
    # print(Δ)
