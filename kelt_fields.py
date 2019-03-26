#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# A module listing the locations of the KELT survey fields.
#
# Rob Siverd
# Created:       2018-05-22
# Last modified: 2019-03-16
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.2.5"

## Python version-agnostic module reloading:
try:
    reload                              # Python 2.7
except NameError:
    try:
        from importlib import reload    # Python 3.4+
    except ImportError:
        from imp import reload          # Python 3.0 - 3.3

## Modules:
import numpy as np
#from numpy.lib.recfunctions import append_fields
#from collections import OrderedDict

## Rectangles on spheres:
import skysquares
reload(skysquares)
sq = skysquares

## FOV rotation module:
try:
    import fov_rotation
    reload(fov_rotation)
except ImportError:
    sys.stderr.write("Module fov_rotation not found!\n")
    sys.exit(1)
rfov = fov_rotation.RotateFOV()

##--------------------------------------------------------------------------##
##------------------         KELT Survey Fields:            ----------------##
##--------------------------------------------------------------------------##

## Each entry consists of (RA_deg, DE_deg, name, active, isdupe)

field_list = []

## Original 13 (North):
field_list.append((  1.50, 31.6656, 'KN01', True, True))
field_list.append(( 30.52, 31.6656, 'KN02', True, True))
field_list.append(( 59.54, 31.6656, 'KN03', True, True))
field_list.append(( 88.56, 31.6656, 'KN04', True, True))
field_list.append((117.58, 31.6656, 'KN05', True, True))
field_list.append((146.60, 31.6656, 'KN06', True, True))
field_list.append((175.62, 31.6656, 'KN07', True, True))
field_list.append((204.64, 31.6656, 'KN08', True, True))
field_list.append((233.66, 31.6656, 'KN09', True, True))
field_list.append((262.68, 31.6656, 'KN10', True, True))
field_list.append((291.70, 31.6656, 'KN11', True, True))
field_list.append((320.72, 31.6656, 'KN12', True, True))
field_list.append((349.74, 31.6656, 'KN13', True, True))

## Joint fields:
field_list.append((114.90,  3.1000, 'KN14', True, True))
field_list.append((253.00,  3.0300, 'KN15', True, True))

## Northern strip:
field_list.append((  0.75, 57.0, 'KN16', True, True))
field_list.append(( 40.75, 57.0, 'KN17', True, True))
field_list.append(( 80.75, 57.0, 'KN18', True, True))
field_list.append((120.75, 57.0, 'KN19', True, True))
field_list.append((160.75, 57.0, 'KN20', True, True))
field_list.append((200.75, 57.0, 'KN21', True, True))
field_list.append((240.75, 57.0, 'KN22', True, True))
field_list.append((280.75, 57.0, 'KN23', True, True))
field_list.append((320.75, 57.0, 'KN24', True, True))

## Polar fields:
field_list.append(( 54.0, 79.0, 'KN25', True, True))
field_list.append((126.0, 79.0, 'KN26', True, True))
field_list.append((198.0, 79.0, 'KN27', True, True))
field_list.append((270.0, 79.0, 'KN28', True, True))
field_list.append((342.0, 79.0, 'KN29', True, True))

## Equatorial stripe:
field_list.append((  6.0, 7.0, 'KN30', True, True))
field_list.append(( 30.0, 7.0, 'KN31', True, True))
field_list.append(( 54.0, 7.0, 'KN32', True, True))
field_list.append(( 78.0, 7.0, 'KN33', True, True))
field_list.append((102.0, 7.0, 'KN34', True, True))
field_list.append((126.0, 7.0, 'KN35', True, True))
field_list.append((150.0, 7.0, 'KN36', True, True))
field_list.append((174.0, 7.0, 'KN37', True, True))
field_list.append((198.0, 7.0, 'KN38', True, True))
field_list.append((222.0, 7.0, 'KN39', True, True))
field_list.append((246.0, 7.0, 'KN40', True, True))
field_list.append((270.0, 7.0, 'KN41', True, True))
field_list.append((294.0, 7.0, 'KN42', True, True))
field_list.append((318.0, 7.0, 'KN43', True, True))
field_list.append((342.0, 7.0, 'KN44', True, True))

##--------------------------------------------------------------------------##
## KELT-South:
field_list.append((  0.00,   3.000, 'KS01', False,  True))
field_list.append(( 23.00,   3.000, 'KS02', False,  True))
field_list.append(( 46.00,   3.000, 'KS03', False,  True))
field_list.append(( 69.00,   3.000, 'KS04', False,  True))

field_list.append(( 91.80,   3.000, 'KS05',  True,  True))
field_list.append((114.90,   3.000, 'KS06',  True,  True)) # GP/2x cadence,
                                                           # joint (KN14)

field_list.append((138.00,   3.000, 'KS07', False,  True))
field_list.append((161.00,   3.000, 'KS08', False,  True))
field_list.append((184.00,   3.000, 'KS09', False,  True))
field_list.append((207.00,   3.000, 'KS10', False,  True))
field_list.append((230.00,   3.000, 'KS11', False,  True))

field_list.append((253.00,   3.000, 'KS12',  True,  True)) # GP/2x cadence,
                                                           # joint (KN15)
field_list.append((276.00,   3.000, 'KS13',  True,  True))
field_list.append((299.00,   3.000, 'KS14',  True,  True))

field_list.append((322.00,   3.000, 'KS15', False,  True))
field_list.append((345.00,   3.000, 'KS16', False,  True))

field_list.append((  0.00, -53.000, 'KS17',  True,  True))
field_list.append(( 23.00, -53.000, 'KS18',  True,  True))
field_list.append(( 46.00, -53.000, 'KS19',  True,  True))
field_list.append(( 69.00, -53.000, 'KS20',  True,  True))
field_list.append(( 91.80, -53.000, 'KS21',  True,  True))
field_list.append((138.00, -20.000, 'KS22',  True,  True))
field_list.append((161.00, -20.000, 'KS23',  True,  True))
field_list.append((184.00, -30.000, 'KS24',  True,  True))
field_list.append((207.00, -30.000, 'KS25',  True,  True))
field_list.append((230.00, -20.000, 'KS26',  True,  True))
field_list.append((299.00, -53.000, 'KS27',  True,  True))
field_list.append((322.00, -53.000, 'KS28',  True,  True))
field_list.append((345.00, -53.000, 'KS29',  True,  True))
field_list.append(( 45.00, -36.000, 'KS30',  True,  True))
field_list.append((  0.85, -29.833, 'KS31', False, False))  # old Blanco1
field_list.append((  0.85, -29.833, 'KS32',  True,  True))  # new Blanco1
field_list.append((123.40, -54.000, 'KS33', False, False))  # old Target1
field_list.append((123.40, -54.000, 'KS34',  True,  True))  # new Target1
field_list.append((115.00, -20.000, 'KS35',  True,  True))
field_list.append((261.00, -53.000, 'KS36',  True,  True))
field_list.append((227.00, -53.000, 'KS37',  True,  True))
field_list.append((192.80, -53.000, 'KS38',  True,  True))
field_list.append((158.64, -53.000, 'KS39',  True,  True))
#field_list.append((261.00, -53.000, 'KS40', kscf_kw))

## KS TEST FIELDS:
field_list.append((255.50, -20.000, 'KS40',  True,  True))
field_list.append((281.00, -20.000, 'KS41',  True,  True))
field_list.append((306.50, -20.000, 'KS42',  True,  True))
field_list.append((332.00, -20.000, 'KS43',  True,  True))
field_list.append(( 90.00, -20.000, 'KS44',  True,  True))

##--------------------------------------------------------------------------##
## Other formats:
_field_ra, _field_de, _field_name, _is_active, _is_dupe = zip(*field_list)

##--------------------------------------------------------------------------##
## Field list as numpy arrays:
def as_arrays():
    return (np.array(_field_ra), np.array(_field_de),
            np.array(_field_name), np.array(_is_active), np.array(_is_dupe))

## Field list as structured array:
def as_structured():
    return np.core.records.fromarrays(as_arrays(), 
                                    names='ra,de,name,active,isdupe')

##--------------------------------------------------------------------------##
## Version reporting:
def get_version():
    return __version__

##--------------------------------------------------------------------------##
## Field list as sky squares:
def borders():
    #kfborders = OrderedDict()
    #kb_radians = [sq.kfield_sky(*coo) for coo in zip(_field_ra, _field_de)]
    diam = 26.2
    kb_radians = [skysquares.kfield_sky(*coo) \
            for coo in zip(_field_ra, _field_de)]
    kb_degrees = [180.0/np.pi*x for x in kb_radians]
    return kb_degrees

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## KELT field point-in-polygon test class:
class FieldChecker(object):

    def __init__(self):
        # hypothetical border for KELT field at RA=0, DE=0:
        self._kf_origin_rad = skysquares.kfield_sky(0.0, 0.0)
        self._kf_origin_deg = np.degrees(self._kf_origin_rad)
        return

    # Compute the difference in orientation (PA) of two vectors. This is the
    # angle subtended by the line segment connecting the "ends" of the two
    # input vectors.
    @staticmethod
    def calc_rel_angle(delta_1, delta_2):
        twopi = 2.0 * np.pi
        dx1, dy1 = delta_1
        dx2, dy2 = delta_2
        theta1 = np.arctan2(dy1, dx1)
        theta2 = np.arctan2(dy2, dx2)
        dtheta = (theta2 - theta1) % twopi
        dtheta[(dtheta > np.pi)] -= twopi
        return dtheta

    # Tally the total angle subtended by border segments. Non-zero totals
    # indicate that the reference point is enclosed by the border. Inputs are:
    # * point --> the RA,DE point to be tested
    # * fcenter --> the RA,DE center of the KELT field to check against
    def count_windings(self, point, fcenter):
        twopi = 2.0 * np.pi
        # de-rotate test point (put test field at origin):
        old_fov = (fcenter[0], fcenter[1], 0.0)
        #ary_pnt = (np.array(point[0]), np.array(point[1]))
        use_point = rfov.roll_to_origin_deg(old_fov, point)
        pra, pde = use_point
        bra, bde = self._kf_origin_deg
        diffs_ra = bra - pra
        diffs_de = bde - pde
        diffs_1 = np.vstack((diffs_ra, diffs_de))
        diffs_2 = np.roll(diffs_1, 1, axis=1)
        windings = np.sum(self.calc_rel_angle(diffs_1, diffs_2)) / twopi
        return np.round(np.abs(windings), 5)


##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Dump border points to file for inspection:
def dump_border(filename, border_pts):
    bra, bde = border_pts
    ra_diff = np.roll(bra, -1) - bra
    de_diff = np.roll(bde, -1) - bde
    sys.stderr.write("bra.shape: %s\n" % str(bra.shape))
    sys.stderr.write("bde.shape: %s\n" % str(bde.shape))
    sys.stderr.write("ra_diff.shape: %s\n" % str(ra_diff.shape))
    sys.stderr.write("de_diff.shape: %s\n" % str(de_diff.shape))
    with open(filename, 'w') as f:
        #for ra,de in border_pts.T:
        for data in zip(bra, bde, ra_diff, de_diff):
            f.write("%12.6f %12.6f %12.6f %12.6f\n" % data)
    return

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##


######################################################################
# CHANGELOG (kelt_fields.py):
#---------------------------------------------------------------------
#
#  2019-03-16:
#     -- Increased __version__ to 0.2.5.
#     -- Added new 'isdupe' column to data set. This is used to indicate
#           if a field is a unique sky position or not. In practice, fields
#           marked isdupe=True should not be plotted. This parameter is
#           needed for Blanco1 and Target1 fields which have separate
#           commissioning/survey designations for identical pointings.
#
#  2019-03-15:
#     -- Increased __version__ to 0.2.0.
#     -- Added as_structured() routine to give field data in structured array.
#
#  2018-09-07:
#     -- Increased __version__ to 0.1.6.
#     -- Added KS40 through KS44 as test fields.
#
#  2018-05-23:
#     -- Increased __version__ to 0.1.5.
#     -- Added FieldChecker() class with point-in-polygon testing capability.
#
#  2018-05-22:
#     -- Increased __version__ to 0.1.0.
#     -- First created kelt_fields.py.
#
