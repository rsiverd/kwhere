#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
#
#    Some useful functions for angle manipulation.
#
# Rob Siverd
# Created:       2011-04-25
# Last modified: 2018-07-19
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "1.7.7"

## Required modules:
import numpy as np
import sys

##--------------------------------------------------------------------------##
## Angle manipulation class:
#class angle:

    ### Initialization:
    #def __init__(self):
    #   pass

    ##-----------------------------------------------------------------------##
    ## Angle reduction functions:
   
# Reduce angle (degrees) to 0 < ang < 360:
def SmallDeg(dAngle):
    """
    Reduce input angle (degrees) to 0 < angle < 360.
    """
    return (360.0 * ((dAngle / 360.0) - np.floor(dAngle / 360.0)))

# Reduce angle (radians) to 0 < ang < 2pi:
def SmallRad(rAngle):
    """
    Reduce input angle (radians) to 0 < angle < 2pi.
    """
    twopi = 2.0 * np.pi
    return (twopi * ((rAngle / twopi) - np.floor(rAngle / twopi)))

# Reduce angle (hours) to 0 < ang < 24:
def SmallHour(hAngle):
    """
    Reduce input angle (hours) to 0 < angle < 24.
    """
    return (24.0 * ((hAngle / 24.0) - np.floor(hAngle / 24.0)))

##-----------------------------------------------------------------------##
## Dimensions-checker:
#def _calc_result_dims(ra1, de1, ra2, de2):
def _coord_dims_okay(ra1, de1, ra2, de2):
    tra1, tde1 = np.atleast_1d(ra1), np.atleast_1d(de1)
    tra2, tde2 = np.atleast_1d(ra2), np.atleast_1d(de2)
    if (tra1.size != tde1.size):
        sys.stderr.write("Input dimension mismatch: ra1, de1\n")
        return False
        #return None
    if (tra2.size != tde2.size):
        sys.stderr.write("Input dimension mismatch: ra2, de2\n")
        return False
        #return None
    if (tra1.size != tra2.size) and (tra1.size > 1) and (tra2.size > 1):
        sys.stderr.write("Incompatible dimensions detected!\n")
        sys.stderr.write("ra1.size == de1.size == %s\n" % tra1.size)
        sys.stderr.write("ra2.size == de2.size == %s\n" % tra2.size)
        return False
        #return None
    #return max(tra1.size, tra2.size)
    return True

## Compute angular separation (radians):
def rAngSep(ra1r, dec1r, ra2r, dec2r): #, safe=True):
    """
    Compute angular separation(s) with a dot product.  All input/output is
    in radians.  Inputs are converted to Cartesian coordinates and their
    dot product is computed.  The arccosine of the dot product is the
    angular separation (since A dot B = |A||B| * cos(angular_separation).
    """

    # Figure out dimensions:
    if not _coord_dims_okay(ra1r, dec1r, ra2r, dec2r):
        return None
    #result_size = _calc_result_dims(ra1r, dec1r, ra2r, dec2r)
    #if result_size == None:
    #    return None

    # Angular differences:
    equal = (ra1r == ra2r) & (dec1r == dec2r)
    #sys.stderr.write("equal: %s\n" % str(equal))
    #sys.stderr.write("equal.dtype: %s\n" % equal.dtype)
    #sys.stderr.write("which: %s\n" % str(which))
    angsep = np.zeros_like(equal, dtype='float')
    #x1 = np.cos(dec1r) * np.cos(ra1r)
    #y1 = np.cos(dec1r) * np.sin(ra1r)
    #z1 = np.sin(dec1r)
    #x2 = np.cos(dec2r) * np.cos(ra2r)
    #y2 = np.cos(dec2r) * np.sin(ra2r)
    #z2 = np.sin(dec2r)
    #dot = x1*x2 + y1*y2 + z1*z2
    dot = np.sin(dec1r) * np.sin(dec2r) \
            + np.cos(dec1r) * np.cos(dec2r) * np.cos(ra1r - ra2r)
    #sys.stderr.write("dot: %s\n" % str(dot))
    #oob = (dot < -1) | (1 < dot)
    #sys.stderr.write("oob: %s\n" % str(oob))
    #angsep[~oob] = np.arccos(dot[~oob])
    #sys.stderr.write("angsep[~equal]: %s\n" % str(angsep[~equal]))
    #sys.stderr.write("dot[~equal]: %s\n" % str(dot[~equal]))
    #angsep[~equal] = np.arccos(dot[~equal])
    #angsep[~equal] = np.arccos(dot[~equal])
    angsep[~equal] = np.arccos(dot[~equal])
    return angsep
    #return np.arccos(dot)
    # ALTERNATIVE:
    # dot = np.sin(dec1r) * np.sin(dec2r) +
    #       np.cos(dec1r) * np.cos(dec2r) * np.cos(ra1r - ra2r)

## Angular separation in degrees (a wrapper for the above):
def dAngSep(ra1d, dec1d, ra2d, dec2d): #, safe=True):
    """
    Compute angular separation(s) using a dot product.  This is a wrapper
    for the rAngSep() function.  See its docstring for more info.
    """
    ra1r, dec1r = np.radians(ra1d), np.radians(dec1d)
    ra2r, dec2r = np.radians(ra2d), np.radians(dec2d)
    return np.degrees(rAngSep(ra1r, dec1r, ra2r, dec2r)) #, safe=safe))

##-----------------------------------------------------------------------##
## Convert Azm/Alt/Lat to HA/Dec (radians, azm reckoned E from N):
def rAzmAltLat_2_HADec(rAzm, rAlt, rLat):
    """
    Convert (azimuth, altitude, latitude) --> (hour angle, declination).
    Inputs must be given in RADIANS.
    Output will be given in RADIANS.
    """
    sys.stderr.write("WARNING: this needs N/S sanity check!\n")
    # Calculate HA:
    numer = np.sin(rAzm)
    denom = np.cos(rAzm)*np.sin(rLat) + np.tan(rAlt)*np.cos(rLat)
    rHA   = np.arctan2(numer, denom)
    # Calculate Dec:
    dummy = np.sin(rLat)*np.sin(rAlt) \
         - np.cos(rLat)*np.cos(rAlt)*np.cos(rAzm)
    rDec  = np.arcsin(dummy)
    return (rHA, rDec)

## Convert Azm/Alt/Lat to HA/Dec (degrees, azm reckoned E from N):
def dAzmAltLat_2_HADec(dAzm, dAlt, dLat):
    """
    Convert (azimuth, altitude, latitude) --> (hour angle, declination).
    Inputs must be given in DEGREES.
    Output will be given in DEGREES.
    """
    #rAzm, rAlt, rLat = np.radians([dAzm, dAlt, dLat])
    rAzm, rAlt, rLat = np.radians(dAzm), np.radians(dAlt), np.radians(dLat)
    #print rAzm,rAlt,rLat
    rHA, rDec = rAzmAltLat_2_HADec(rAzm, rAlt, rLat)
    return (np.degrees(rHA), np.degrees(rDec))
    #return np.degrees(rAzmAltLat_2_HADec(rAzm, rAlt, rLat))

##-----------------------------------------------------------------------##
##-----------------------------------------------------------------------##

##-----------------------------------------------------------------------##
## Convert HA/Dec/Lat to Azm/Alt (radians, Azm *WEST from SOUTH*):
def rHADecLat_2_sAzAlt(rHA, rDec, rLat):
    """
    Converts (hour angle, declination, latitude) --> (azimuth, altitude).
    -- Azimuth is reckoned *WEST from SOUTH*.
    -- Latitude is measured positively northwards.
    -- Hour Angle is measured positively westwards.
    -- Inputs must be given in RADIANS.
    -- Output will be given in RADIANS.
    """
    # Calculate Azm:
    numer = np.sin(rHA)
    denom = np.cos(rHA)*np.sin(rLat) - np.tan(rDec)*np.cos(rLat)
    rAzm = np.arctan2(numer, denom)
    # Calculate Alt:
    dummy = np.sin(rLat)*np.sin(rDec) \
         + np.cos(rLat)*np.cos(rDec)*np.cos(rHA)
    rAlt  = np.arcsin(dummy)
    return (rAzm, rAlt)

##-----------------------------------------------------------------------##
## Convert HA/Dec/Lat to Azm/Alt (degrees, Azm WEST from SOUTH):
def dHADecLat_2_sAzAlt(dHA, dDec, dLat):
    """
    Converts (hour angle, declination, latitude) --> (azimuth, altitude).
    -- Azimuth is reckoned *WEST from SOUTH*.
    -- Latitude is measured positively northwards.
    -- Hour Angle is measured positively westwards.
    -- Inputs must be given in DEGREES.
    -- Output will be given in DEGREES.
    """
    rHA, rDec, rLat = np.radians(dHA), np.radians(dDec), np.radians(dLat)
    rAz, rAlt = rHADecLat_2_sAzAlt(rHA, rDec, rLat)
    return (np.degrees(rAz), np.degrees(rAlt))

##-----------------------------------------------------------------------##
## Convert HA/Dec/Lat to Azm/Alt (radians, Azm EAST from NORTH):
def rHADecLat_2_nAzAlt(rHA, rDec, rLat):
    """
    Converts (hour angle, declination, latitude) --> (azimuth, altitude).
    -- Azimuth is reckoned *EAST from NORTH*.
    -- Latitude is measured positively northwards.
    -- Hour Angle is measured positively westwards.
    -- Inputs must be given in RADIANS.
    -- Output will be given in RADIANS.
    """
    rAz, rAlt = rHADecLat_2_sAzAlt(rHA, rDec, rLat)
    rAz = (rAz + np.pi) % np.radians(360.0)
    return (rAz, rAlt)

##-----------------------------------------------------------------------##
## Convert HA/Dec/Lat to Azm/Alt (degrees, Azm EAST from NORTH):
def dHADecLat_2_nAzAlt(dHA, dDec, dLat):
    """
    Converts (hour angle, declination, latitude) --> (azimuth, altitude).
    -- Azimuth is reckoned *EAST from NORTH*.
    -- Latitude is measured positively northwards.
    -- Hour Angle is measured positively westwards.
    -- Inputs must be given in DEGREES.
    -- Output will be given in DEGREES.
    """
    # Unit conversion:
    rHA, rDec, rLat = np.radians(dHA), np.radians(dDec), np.radians(dLat)
    rAz, rAlt = rHADecLat_2_nAzAlt(rHA, rDec, rLat)
    return (np.degrees(rAz), np.degrees(rAlt))

##--------------------------------------------------------------------------##
##*********************  Spherical Location Estimates: *********************##
##--------------------------------------------------------------------------##

## Average direction of vectors on unit sphere (RADIANS):
def spheremean_rad(RA_rad, DE_rad, dev=False):
    """
    Compute mean (RA, Dec) for a set of RA, Dec directions (RADIANS).
    Returns:
         (avg_RA, avg_DE)               # dev=False
         (avg_RA, avg_DE, angvar)       # dev=True
    """
    vecX = np.sum(np.cos(DE_rad) * np.cos(RA_rad))     # total X length
    vecY = np.sum(np.cos(DE_rad) * np.sin(RA_rad))     # total Y length
    vecZ = np.sum(np.sin(DE_rad))                      # total Z length
    R_tot = np.sqrt(vecX*vecX + vecY*vecY + vecZ*vecZ) # total distance
    angvar = 1.0 - R_tot                               # 'circular' variance
    avg_DE = np.arcsin(vecZ / R_tot)
    avg_RA = np.arctan2(vecY, vecX) % (2.0 * np.pi)
    if dev:
        return (avg_RA, avg_DE, angvar)
    else:
        return (avg_RA, avg_DE)

## Average direction of vectors on unit sphere (RADIANS):
def spheremean_deg(RA_deg, DE_deg, dev=False):
    """
    Compute spherical mean and angular variance for a set of RA, Dec 
    coordinates (DEGREES).
 
    Returns:
         (avg_RA, avg_DE)               # dev=False
         (avg_RA, avg_DE, angvar)       # dev=True
    """
    RA_rad = np.radians(RA_deg)
    DE_rad = np.radians(DE_deg)
    avg_RA_r, avg_DE_r, angvar = spheremean_rad(RA_rad, DE_rad, dev=True)
    if dev:
        return (np.degrees(avg_RA_r), np.degrees(avg_DE_r), angvar)
    else:
        return (np.degrees(avg_RA_r), np.degrees(avg_DE_r))

## Medoid direction and scatter of vectors on unit sphere (RADIANS):
def sphere_medoid_rad(RA_rad, DE_rad, dev=False):
    """
    Find medoid of a set of RA, Dec coordinates (RADIANS). Optionally returns
    sum of absolute differences of input coordinates from medoid.

    Returns:
         (med_RA, med_DE)                   # dev=False
         (med_RA, med_DE, sum_abs_diffs)    # dev=True
    """
    sad = np.zeros_like(RA_rad)  # sum of absolute diffs
    for i, (try_RA, try_DE) in enumerate(zip(RA_rad, DE_rad)):
       sad[i] = np.sum(rAngSep(RA_rad, DE_rad, try_RA, try_DE))
    mid = sad.argmin()
    if dev:
        return (RA_rad[mid], DE_rad[mid], sad[mid])
    else:
        return (RA_rad[mid], DE_rad[mid])

## Medoid direction and scatter of vectors on unit sphere (DEGREES):
def sphere_medoid_deg(RA_deg, DE_deg, dev=False):
    """
    Find medoid of a set of RA, Dec coordinates (DEGREES). Optionally returns
    sum of absolute differences of input coordinates from medoid.

    Returns:
         (med_RA, med_DE)                   # dev=False
         (med_RA, med_DE, sum_abs_diffs)    # dev=True
    """
    RA_rad = np.radians(RA_deg)
    DE_rad = np.radians(DE_deg)
    med_RA_r, med_DE_r, sad = sphere_medoid_rad(RA_rad, DE_rad, dev=True)
    if dev:
        return (np.degrees(med_RA_r), np.degrees(med_DE_r), sad)
    else:
        return (np.degrees(med_RA_r), np.degrees(med_DE_r))



######################################################################
# CHANGELOG (angle.py):
#---------------------------------------------------------------------
#
#  2018-07-19:
#     -- Increased __version__ to 1.7.7.
#     -- Degree form of AzAltLat->HADec now returns a tuple like its radian
#           counterpart.
#     -- Simplified radian conversion in AzAltLat->HADec converter.
#
#  2018-02-12:
#     -- Increased __version__ to 1.7.6.
#     -- Fixed check for bad/inconsistent dimensions used in AngSep routines.
#           Previously, the (allowed) case of all dimensions equal was not
#           counted as valid. Now it should work.
#
#  2018-02-03:
#     -- Increased __version__ to 1.7.5.
#     -- Results of rAngSep are now initialized to size of 'equal' array with
#           explicit float type. This guarantees the proper type/dimensions.
#     -- Added dimensionality check routine to detect bogus rAngSep inputs.
#
#  2018-01-25:
#     -- Increased __version__ to 1.7.0.
#     -- rAngSep and dAngSep now explicitly handle cases of equal inputs to
#           ensure 0-valued results and avoid NaNs in response.
#     -- Added sphere_medoid_deg() method.
#     -- Indentation is now 4 spaces.
#
#  2017-12-13:
#     -- Increased __version__ to 1.6.5.
#     -- Now convert each input to radians separately in dAngSep(). This
#           ensures that dimensionality and axes are preserved (the old
#           syntax did not seem to do this with multi-dimensional arrays).
#
#  2016-11-08:
#     -- Increased __version__ to 1.6.0.
#     -- Fixed array support in dHADecLat_2_nAzAlt().
#
#  2015-01-02:
#     -- Increased __version__ to 1.5.0.
#     -- Added sphere_medoid_rad() 
#     -- Added spheremean_rad() and spheremean_deg() directional averages.
#
#  2014-04-28:
#     -- Increased __version__ to 1.0.1.
#     -- Now explicitly specify UTF-8 file encoding.
#
#  2012-09-26:
#     -- Added docstrings to several functions.
#
#  2012-06-27:
#     -- Added azm/alt/lat --> HA/dec conversion functions.
#     -- Added HA/dec/lat --> azm/alt conversion functions.
#     -- Improved docstring formatting.
#
#  2011-04-25:
#     -- Added SmallDeg, SmallRad, & SmallHour functions.
#     -- First created angle.py.
#
