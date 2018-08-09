#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
#
#    Test codes to project squares onto celestial sphere, i.e., given the
# center and dimensions of a rectangular FOV, directly compute RA,Dec border
# coordinates.
#
# Rob Siverd
# Created:       2013-12-13
# Last modified: 2018-05-23
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "1.0.5"

## Modules:
import os
import sys
import time
import numpy as np

##--------------------------------------------------------------------------##

def rade2xyz(rad_ra, rad_de):
    x = np.cos(rad_de) * np.cos(rad_ra)
    y = np.cos(rad_de) * np.sin(rad_ra)
    z = np.sin(rad_de)
    return np.vstack((x, y, z))

def deg2xyz(deg_ra, deg_de):
    return rade2xyz(np.radians(deg_ra), np.radians(deg_de))

## Convert RA, Dec to Cartesian:
def equ2xyz(equ_pts):
    ra = equ_pts[0]
    de = equ_pts[1]
    x = np.cos(de) * np.cos(ra)
    y = np.cos(de) * np.sin(ra)
    z = np.sin(de)
    return np.vstack((x, y, z))

## Convert Cartesian points to RA, Dec:
def xyz2equ(xyz_pts):
    # Shape/dimension sanity check:
    if ((xyz_pts.ndim != 2) or (xyz_pts.shape[0] != 3)):
        sys.stderr.write("XYZ points have wrong shape!\n")
        return (0,0)
    tx = np.array(xyz_pts[0]).flatten()
    ty = np.array(xyz_pts[1]).flatten()
    tz = np.array(xyz_pts[2]).flatten()
    ra = np.arctan2(ty, tx)
    de = np.arcsin(tz)
    equ_coo = np.vstack((ra, de))
    return equ_coo

## FIXME FIXME FIXME --> switch to skyrotation module   FIXME FIXME FIXME
##--------------------------------------------------------------------------##
## Rotation matrices:
def calc_xrot_rad(ang):
    return np.matrix([[  1.0,          0.0,          0.0],
                      [  0.0,  np.cos(ang), -np.sin(ang)],
                      [  0.0,  np.sin(ang),  np.cos(ang)]])

def calc_yrot_rad(ang):
    return np.matrix([[ np.cos(ang),   0.0,  np.sin(ang)],
                      [         0.0,   1.0,          0.0],
                      [-np.sin(ang),   0.0,  np.cos(ang)]])

def calc_zrot_rad(ang):
    return np.matrix([[ np.cos(ang), -np.sin(ang),   0.0],
                      [ np.sin(ang),  np.cos(ang),   0.0],
                      [         0.0,          0.0,   1.0]])

## Degree versions:
def calc_xrot_deg(angle):
    return calc_xrot_rad(np.radians(angle))

def calc_yrot_deg(angle):
    return calc_yrot_rad(np.radians(angle))

def calc_zrot_deg(angle):
    return calc_zrot_rad(np.radians(angle))

##--------------------------------------------------------------------------##
## How to rotate an array of vectors:
def rotate_xyz(rmatrix, xyz_list):
    #result = np.array([np.dot(i, rmatrix.T) for i in xyz_list])
    #result = np.array([np.dot(i, rmatrix.T) for i in xyz_list])
    #return np.squeeze(result)
    return np.dot(rmatrix, xyz_list)

def xrotate_xyz(ang_deg, xyz_list):
    rot_mat = calc_xrot_deg(ang_deg)
    return np.dot(rot_mat, xyz_list)

def yrotate_xyz(ang_deg, xyz_list):
    rot_mat = calc_yrot_deg(ang_deg)
    return np.dot(rot_mat, xyz_list)

def zrotate_xyz(ang_deg, xyz_list):
    rot_mat = calc_zrot_deg(ang_deg)
    return np.dot(rot_mat, xyz_list)

def xrot(ang_deg, xyz_list):
    return xrotate_xyz(ang_deg, xyz_list)

def yrot(ang_deg, xyz_list):
    return yrotate_xyz(ang_deg, xyz_list)

def zrot(ang_deg, xyz_list):
    return zrotate_xyz(ang_deg, xyz_list)

##--------------------------------------------------------------------------##
## Shift *great cirlce* points into (-pi < ang < pi):
def adjust_GC_RA(coo_vec):
    # fix RA range:
    ra_pts = coo_vec[0]
    ra_pts = ra_pts % (2.0 * np.pi)
    ra_pts -= np.pi
    # sort coords by RA:
    return coo_vec[:, coo_vec[0].argsort()]

##--------------------------------------------------------------------------##
## Define equatorial great circle (equator):
def make_egc(pts=100, arc=1.0):
    #t_ra = arc * np.pi * np.linspace(-1.0, 1.0, pts)
    t_ra = arc * np.pi * np.linspace(-1.0, 1.0, num=pts, endpoint=True)
    t_de = np.zeros_like(t_ra)
    return np.vstack((t_ra, t_de))

def test_make_egc(pts=100, arc_deg=360.0):
    #t_ra = arc * np.pi * np.linspace(-1.0, 1.0, pts)
    t_ra = np.radians(arc_deg * np.linspace(-0.5, 0.5, num=pts, endpoint=True))
    t_de = np.zeros_like(t_ra)
    return np.vstack((t_ra, t_de))

# --------
def make_mer(pts=100, arc=1.0):
    t_de = arc * np.pi * np.linspace(-0.5, 0.5, num=pts, endpoint=True)
    t_ra = np.zeros_like(t_de)
    return np.vstack((t_ra, t_de))

def test_make_mer(pts=100, arc_deg=180.0, ndrop=2):
    #t_de = arc * np.linspace(-90.0, 90.0, pts)
    t_de = np.radians(arc_deg * np.linspace(-0.5, 0.5, num=pts, endpoint=True))
    t_ra = np.zeros_like(t_de)
    return np.vstack((t_ra, t_de))

## Cartesian equatorial great circle generator:
def cart_egc(**kwargs):
    return equ2xyz(make_egc(**kwargs))

def test_cart_egc(**kwargs):
    return equ2xyz(test_make_egc(**kwargs))

## Cartesian meridian great circle generator:
def cart_mer(**kwargs):
    return equ2xyz(make_mer(**kwargs))

def test_cart_mer(**kwargs):
    return equ2xyz(test_make_mer(**kwargs))

##--------------------------------------------------------------------------##
def merid_pts_rad(rad_ra, pts=50):
    de_vec = np.pi * np.linspace(-0.5, 0.5, pts)
    ra_vec = 0.0 * de_vec + rad_ra
    return np.vstack((ra_vec, de_vec))

#def gcdraw(ax, xyz_pts, c='b', **kwargs):
def gcdraw(ax, xyz_pts, c='b', s=0.1, **kwargs):
    equ_pts = xyz2equ(xyz_pts)
    #ax.plot(equ_pts[0], equ_pts[1], color=color, **kwargs)
    ax.scatter(equ_pts[0], equ_pts[1], color=c, s=s, **kwargs)


## FIXME: this routine needs work ... borders overextend slightly
def new_kborder(pts=150, frac=0.07):
    #frac = 0.07
    #frac = 0.17
    ### 1 extra point ( pts=150) at diam=10.0
    ### 2 extra points (pts=150) at diam=20.0
    diam = 26.2
    #diam =  30.0
    #arcfix = 1.0 - (3.0 / float(pts))
    #arcfix = 1.0
    #plane = test_cart_egc(pts=pts, arc_deg=diam*arcfix)
    #merid = test_cart_mer(pts=pts, arc_deg=diam*arcfix)
    frac = 0.071
    plane = cart_egc(pts=150, arc=1*frac)
    merid = cart_mer(pts=150, arc=2*frac)

    #ft_xyz = yrot(-13.0, xrot(  0.0, plane))  # top
    #fl_xyz = zrot( 13.0, xrot(180.0, merid))  # left
    #fb_xyz = yrot( 13.0, xrot(180.0, plane))  # bottom
    #fr_xyz = zrot(-13.0, xrot(  0.0, merid))  # right
    ft_xyz = yrot(-0.5*diam, xrot(  0.0, plane))  # top
    fl_xyz = zrot( 0.5*diam, xrot(180.0, merid))  # left
    fb_xyz = yrot( 0.5*diam, xrot(180.0, plane))  # bottom
    fr_xyz = zrot(-0.5*diam, xrot(  0.0, merid))  # right
    kf_pts = np.hstack((ft_xyz, fl_xyz, fb_xyz, fr_xyz))

    #ft_xyz = yrot(-13.0, plane)[::+1] # top
    #fl_xyz = zrot( 13.0, merid)[::+1] # left
    #fb_xyz = yrot( 13.0, plane)[::-1] # bottom
    #fr_xyz = zrot(-13.0, merid)[::-1] # right
    #kf_pts = np.hstack((ft_xyz, fl_xyz, fb_xyz, fr_xyz))
    return kf_pts

## Generate KELT field border points (Cartesian):
def kfield_xyz(ra_deg, de_deg):
    diam = 26.2
    #kb_pts = new_skyrect(diam, diam)
    kb_pts = new_kborder()
    kb_pts = yrot(-de_deg, kb_pts) # rotate in Dec
    kb_pts = zrot( ra_deg, kb_pts) # rotate in RA
    return kb_pts

## Generate KELT field border points (RA, DE):
def kfield_sky(ra_deg, de_deg):
    return xyz2equ(kfield_xyz(ra_deg, de_deg))

##--------------------------------------------------------------------------##
## Viewing region for basemap (lon/lat corners):
def rect_corners_xyz(diam_deg):
    plane = test_cart_egc(pts=2, arc_deg=diam_deg)
    top   = yrot(-0.5*diam_deg, xrot(  0.0, plane))
    bot   = yrot( 0.5*diam_deg, xrot(180.0, plane))
    return np.hstack((top, bot))

def lonlat_corners_rad(diam_deg):
    return xyz2equ(rect_corners_xyz(diam_deg))

def lonlat_corners_deg(diam_deg):
    return np.degrees(lonlat_corners_rad(diam_deg))

##--------------------------------------------------------------------------##
def derp_skysquare(diam_deg, pts=150):
    plane = test_cart_egc(pts=pts, arc_deg=diam_deg)
    upper = yrot(-0.5*diam_deg, xrot(180.0, plane))     # "top" of square
    rlist = 90.0 * np.arange(4)
    sides = [xrot(rr, upper) for rr in rlist]
    return xyz2equ(np.hstack(sides))

def new_skyrect(width, height, pts=150):
    wfrac = width / 360.0
    hfrac = height / 180.0
    plane = cart_egc(pts=pts, arc=wfrac)
    merid = cart_mer(pts=pts, arc=hfrac)
    ft_xyz = yrot(-0.5*height, xrot(  0.0, plane))  # top
    fl_xyz = zrot( 0.5*width,  xrot(180.0, merid))  # left
    fb_xyz = yrot( 0.5*height, xrot(180.0, plane))  # bottom
    fr_xyz = zrot(-0.5*width,  xrot(  0.0, merid))  # right
    return np.hstack((ft_xyz, fl_xyz, fb_xyz, fr_xyz))

def test_skyrect_xyz(width_deg, height_deg, pts=150):
    plane = test_cart_egc(pts=pts, arc_deg=width_deg)
    merid = test_cart_mer(pts=pts, arc_deg=height_deg)
    sys.stderr.write("plane.shape: %s\n" % str(plane.shape))
    sys.stderr.write("plane:       %s\n" % str(plane))
    sys.stderr.write("merid.shape: %s\n" % str(merid.shape))
    sys.stderr.write("merid:       %s\n" % str(merid))
    ft_xyz = yrot(-0.5*height_deg, xrot(  0.0, plane))  # top
    fl_xyz = zrot( 0.5*width_deg,  xrot(180.0, merid))  # left
    fb_xyz = yrot( 0.5*height_deg, xrot(180.0, plane))  # bottom
    fr_xyz = zrot(-0.5*width_deg,  xrot(  0.0, merid))  # right
    return np.hstack((ft_xyz, fl_xyz, fb_xyz, fr_xyz))

def test_skyrect(width_deg, height_deg, pts=150):
    border_xyz = test_skyrect_xyz(width_deg, height_deg, pts)
    return xyz2equ(border_xyz)

## Generate KELT field border points (RA, DE):
def ssquare_sky(ra_deg, de_deg):
    return xyz2equ(kfield_xyz(ra_deg, de_deg))

## FIXME (this should be removed):
#def new_moon_grid(pts=15, diam_deg=50.0, npts=10):
#   #sys.stderr.write("diam_deg (new_moon_grid): %.3f\n" % diam_deg)
#   plane = cart_egc(pts=npts, arc=1.0*diam_deg / 360.)
#   merid = cart_mer(pts=npts, arc=2.0*diam_deg / 360.)
#   stripes = []
#   for ddec in np.linspace(-0.5, 0.5, npts):
#      #ddec = 0.5 * diam * snum
#      next_row = yrot(diam_deg*ddec, plane)
#      #next_row = yrot(0.5*diam*snum, plane)
#      stripes.append(next_row)
#      #mg_pts = np.hstack((mg_pts, next_row))
#   mg_pts = np.hstack((stripes))
#
#   return mg_pts
#
#def moongrid_xyz(ra_deg, de_deg, diam_deg, pts):
#   #sys.stderr.write("diam_deg (moongrid_xyz): %.3f\n" % diam_deg)
#   #mg_pts = new_moon_grid(diam_deg / 360.0)
#   mg_pts = new_moon_grid(pts, diam_deg)
#   mg_pts = yrot(-de_deg, mg_pts) # rotate in Dec
#   mg_pts = zrot( ra_deg, mg_pts) # rotate in RA
#   return mg_pts

### Generate KELT field border points (RA, DE):
#def moongrid_sky(ra_deg, de_deg, diam_deg, pts=10):
#   return xyz2equ(moongrid_xyz(ra_deg, de_deg, diam_deg, pts))

##--------------------------------------------------------------------------##



##--------------------------------------------------------------------------##
#fig_dims = (16, 10)
#fig = plt.figure(1, figsize=fig_dims)
#fig.clf()
##fig.subplots_adjust(left=0.07, right=0.95)
#ax1 = fig.add_subplot(111, projection='sky_hammer_180')
#ax1.grid(True)



##--------------------------------------------------------------------------##
#mdraw_deg(ax1,  13.0)
#mdraw_deg(ax1, -13.0)

#ax1.plot(eq_ra, eq_de, c='r')

#gcdraw(ax1, gcpts)
#gcdraw(ax1, plane)

#ax1.plot(tra, tde, c='g', lw=10)

#yrdraw(ax1,  13.1, 'g')
#yrdraw(ax1, -13.1, 'g')

#gcdraw(ax1, kfield_xyz( 175.0,  45.0), c='m')
#
#gcdraw(ax1, kfield_xyz( 180.0,   3.0), c='y')
#
#gcdraw(ax1, kfield_xyz( 185.0, -53.0), c='c')
#
#gcdraw(ax1, kfield_xyz(  74.0,  82.0), c='g')
#gcdraw(ax1, kfield_xyz( -66.0,  22.0), c='r')
#gcdraw(ax1, kfield_xyz(  -5.0, -44.0), c='b')
#gcdraw(ax1, kfield_xyz(  74.0,  20.0), c='c')
#gcdraw(ax1, kfield_xyz( -54.0, -85.0), c='orange')
#gcdraw(ax1, kfield_xyz( -90.0, -73.0), c='lime')
#gcdraw(ax1, kfield_xyz(  60.0, -55.0), c='r')

##--------------------------------------------------------------------------##


#plt.tight_layout() # adjust boundaries sensibly, matplotlib v1.1+
#plt.draw()






######################################################################
# CHANGELOG (skysquares.py):
#---------------------------------------------------------------------
#
#  2018-05-23:
#     -- Increased __version__ to 1.0.5.
#     -- Border points are now in provided in self-consistent order from
#           new_kborder() and new_skyrect().
#
#  2018-02-21:
#     -- Increased __version__ to 1.0.1.
#     -- Changed indentation to 4 spaces.
#
#  2014-06-13:
#     -- Increased __version__ to 1.0.0.
#     -- Moved code to skysquares.py module.
#
#  2013-12-13:
#     -- First created square_test.py.
#
