#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Check which field(s) contain the specified RA,DE point on the sky.
#
# Rob Siverd
# Created:       2018-05-22
# Last modified: 2018-09-25
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.2.2"

## Python version-agnostic module reloading:
try:
    reload                              # Python 2.7
except NameError:
    try:
        from importlib import reload    # Python 3.4+
    except ImportError:
        from imp import reload          # Python 3.0 - 3.3

## Modules:
import argparse
import warning
import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#import matplotlib.ticker as mt
#import matplotlib._pylab_helpers as hlp
#from matplotlib.colors import LogNorm
#from matplotlib import colors
#import matplotlib.colors as mplcolors
#import matplotlib.gridspec as gridspec

import angle
reload(angle)

import skysquares
reload(skysquares)
sq = skysquares

## Use basemap for plots, if available:
_have_basemap = False
try:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        from mpl_toolkits.basemap import Basemap
    _have_basemap = True
except ImportError:
    sys.stderr.write("Basemap not installed ... plots will be ugly!\n")


##--------------------------------------------------------------------------##
## KELT fields:
try:
    import kelt_fields
    reload(kelt_fields)
except ImportError:
    sys.stderr.write("Module kelt_fields not found!\n")
    sys.exit(1)
kfc = kelt_fields.FieldChecker()

## FOV rotation module:
try:
    import fov_rotation
    reload(fov_rotation)
except ImportError:
    sys.stderr.write("Module fov_rotation not found!\n")
    sys.exit(1)
rfov = fov_rotation.RotateFOV()

##--------------------------------------------------------------------------##
## Colors for fancy terminal output:
NRED    = '\033[0;31m'   ;  BRED    = '\033[1;31m'
NGREEN  = '\033[0;32m'   ;  BGREEN  = '\033[1;32m'
NYELLOW = '\033[0;33m'   ;  BYELLOW = '\033[1;33m'
NBLUE   = '\033[0;34m'   ;  BBLUE   = '\033[1;34m'
NMAG    = '\033[0;35m'   ;  BMAG    = '\033[1;35m'
NCYAN   = '\033[0;36m'   ;  BCYAN   = '\033[1;36m'
NWHITE  = '\033[0;37m'   ;  BWHITE  = '\033[1;37m'
ENDC    = '\033[0m'

## Suppress colors in cron jobs:
if (os.getenv('FUNCDEF') == '--nocolors'):
    NRED    = ''   ;  BRED    = ''
    NGREEN  = ''   ;  BGREEN  = ''
    NYELLOW = ''   ;  BYELLOW = ''
    NBLUE   = ''   ;  BBLUE   = ''
    NMAG    = ''   ;  BMAG    = ''
    NCYAN   = ''   ;  BCYAN   = ''
    NWHITE  = ''   ;  BWHITE  = ''
    ENDC    = ''

## Fancy text:
degree_sign = u'\N{DEGREE SIGN}'

## Dividers:
halfdiv = "----------------------------------------"
fulldiv = halfdiv + halfdiv

##--------------------------------------------------------------------------##
def ldmap(things):
    return dict(zip(things, range(len(things))))

def argnear(vec, val):
    return (np.abs(vec - val)).argmin()

## Settings:
debug = False
timer = False
vlevel = 0
#prog_name = 'kwhere.py'
prog_name = os.path.basename(__file__)
full_prog = sys.argv[0]
base_prog = os.path.basename(full_prog)
num_todo = 0
plot_radius_deg = 50.0      # field center separation cutoff for plot inclusion
plot_view_diam  = 45.0

##--------------------------------------------------------------------------##
##*********************     Help and options menu:     *********************##
##--------------------------------------------------------------------------##


##--------------------------------------------------------------------------##
## Parse arguments and run script:
if __name__ == '__main__':

    # ------------------------------------------------------------------
    descr_txt = """
    Check which field(s) contain the specified RA/DE point.
    
    Version: %s
    """ % __version__
    parser = argparse.ArgumentParser(
            prog=os.path.basename(__file__),
            description=descr_txt)
    parser.add_argument('ra', help='test coordinate RA', type=float)
    parser.add_argument('de', help='test coordinate DE', type=float)
    parser.add_argument('-P', '--plot', dest='plot_name', default=None,
            help='Produce plot showing point placement near fields')
    parser.add_argument('--debug', dest='debug', default=False,
             help='Enable extra debugging messages', action='store_true')
    parser.add_argument('-q', '--quiet', action='count', default=0)
    parser.add_argument('-v', '--verbose', action='count', default=0)
    #parser.add_argument('remainder', help='other stuff', nargs='*')
    # ------------------------------------------------------------------

    context = parser.parse_args()
    #context.vlevel = context.verbose - context.quiet
    context.vlevel = 99 if context.debug else (context.verbose-context.quiet)

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
## Load KELT field locations and borders:
kra, kde, kname, active = kelt_fields.as_arrays()
kf_borders = kelt_fields.borders()

##--------------------------------------------------------------------------##
## Brute-force borders and winding number calculation:
tik = time.time()
nfields = len(kf_borders)
kfb_sep_stats = {'min':[], 'max':[], 'avg':[], 'wind':[],}
for i,kfedge in enumerate(kf_borders):
    #sys.stderr.write("\rField %d of %d ... " % (i+1, nfields))
    kfbra, kfbde = kfedge
    kfb_sep = angle.dAngSep(kfbra, kfbde, context.ra, context.de)
    fcenter = (kra[i], kde[i])
    winding = kfc.count_windings((context.ra, context.de), fcenter)
    kfb_sep_stats['min'].append(kfb_sep.min())
    kfb_sep_stats['max'].append(kfb_sep.max())
    kfb_sep_stats['avg'].append(kfb_sep.mean())
    kfb_sep_stats['wind'].append(winding)
    pass

tok = time.time()
if (context.vlevel >= 1):
    sys.stderr.write("done. Took %.3f seconds.\n" % (tok-tik))
#n_windings = np.array(n_windings)

## Promote results to numpy arrays:
for kk in kfb_sep_stats.keys():
    kfb_sep_stats[kk] = np.array(kfb_sep_stats[kk])


#windings = np.array(zip(*kfb_sep_stats)[3])
foundidx = kfb_sep_stats['wind'].nonzero()[0]

n_containing = foundidx.size
if (n_containing == 0):
    sys.stdout.write("Test point is not contained by any KELT field.\n")
    sys.exit(0)

sys.stdout.write("Found %d KELT field(s) containing the test point:\n"
        % n_containing)
for ii in foundidx:
    fra, fde, fname = kra[ii], kde[ii], kname[ii]
    sys.stderr.write("%s --- (%7.3f, %+7.3f)\n" % (fname, fra, fde))

## Order by decreasing separation:
n_kept = np.sum(kfb_sep_stats['avg'] <= plot_radius_deg)
fsep_order = np.argsort(kfb_sep_stats['avg'])

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
## Construct viewing region centered around test point:
corners = sq.lonlat_corners_deg(plot_view_diam)        # at RA,DE = 0,0
targfov = (context.ra, context.de, 0.0)
viewbox = rfov.migrate_fov_deg((0., 0., 0.), targfov, corners)

## Latitude lines (parallels):
n_lat_lines = 10
vb_lat_range = viewbox[1].max() - viewbox[1].min()
approx_lat_spacing = 2.0**np.round(np.log2(vb_lat_range / n_lat_lines))
par_list = np.arange( -90.0,  90.0, approx_lat_spacing) # latitude parallels

## Longitude lines (meridians):
n_lon_lines = 5
#vb_lon_range = viewbox[0].max() - viewbox[0].min()
vb_lon_range = plot_view_diam / np.cos(np.radians(context.de))
approx_lon_spacing = 2.0**np.round(np.log2(vb_lon_range / n_lon_lines))
#if approx_lon_spacing > 5:
#if (context.de <= -70.0) or (70.0 <= context.de):
#    approx_lon_spacing = 20.0
if approx_lon_spacing > 20.0:
    approx_lon_spacing = 20.0
mer_list = np.arange(-180.0, 180.0, approx_lon_spacing) # longitude meridians

## Basemap config:
#proj_type = 'hammer'
#proj_type = 'gnom'

## Viewing box corners:
#celestial = False
#if celestial:
#    ll_corners = {
#        'llcrnrlon':viewbox[0][2],
#        'llcrnrlat':viewbox[1][2],
#        'urcrnrlon':viewbox[0][0],
#        'urcrnrlat':viewbox[1][0],}
#else:
#    ll_corners = {
#        'llcrnrlon':viewbox[0][3],
#        'llcrnrlat':viewbox[1][3],
#        'urcrnrlon':viewbox[0][1],
#        'urcrnrlat':viewbox[1][1],}


#proj_kw = {'projection':'hammer', 'lon_0':context.ra, }
proj_kw = {
    'projection':'gnom',
    'lon_0':context.ra,
    'lat_0':context.de,
    'celestial':False,
    'llcrnrlon':viewbox[0][3],
    'llcrnrlat':viewbox[1][3],
    'urcrnrlon':viewbox[0][1],
    'urcrnrlat':viewbox[1][1],
    }

#proj_kw['celestial'] = True
#lon_center = context.ra

##--------------------------------------------------------------------------##
## Plot config:

# gridspec examples:
# https://matplotlib.org/users/gridspec.html

#gs1 = gridspec.GridSpec(4, 4)
#gs1.update(wspace=0.025, hspace=0.05)  # set axis spacing

##--------------------------------------------------------------------------##
fig_dims = (10, 9)
fig = plt.figure(1, figsize=fig_dims)
plt.gcf().clf()
#fig, axs = plt.subplots(2, 2, sharex=True, figsize=fig_dims, num=1)
# sharex='col' | sharex='row'
#fig.frameon = False # disable figure frame drawing
#fig.subplots_adjust(left=0.07, right=0.95)
#ax1 = plt.subplot(gs[0, 0])
ax1 = fig.add_subplot(111)
#ax1 = fig.add_axes([0.07, 0.07, 0.95, 0.95])

## Initialize the sky grid as needed:
if _have_basemap:
    sky = Basemap(ax=ax1, **proj_kw)
    grid_kw = {'linewidth':0.5, 'labelstyle':"+/-"}
    par_labs = [1, 1, 0, 0]     #  latitudes only on left/right
    mer_labs = [0, 0, 1, 1]     # longitudes only on left/right
    sky.drawparallels(par_list, labels=par_labs, **grid_kw)
    sky.drawmeridians(mer_list, labels=mer_labs, **grid_kw)
else:
    sky = None
    ax1.grid(True)

## Show borders of fields within the viewing area:
for idx in fsep_order[:n_kept]:
    kfra, kfde = kf_borders[idx]
    fname = kname[idx]
    if sky:
        fx, fy = sky(kfra, kfde)
        #sky.scatter(fx, fy, lw=0, s=1, label=fname)
        #sky.plot(fx, fy, marker='.', ms=2, label=fname)
        sky.plot(fx, fy, label=fname)
    else:
        ax1.plot(kfra, kfde, marker='.', ls='', ms=1, label=fname)

if sky:
    fx, fy = sky(context.ra, context.de)
    sky.plot(fx, fy, c='r', marker='x', ms=20, ls='')
else:
    ax1.plot(context.ra, context.de, c='r', marker='x', ms=20, ls='')
ax1.legend(loc='upper right')

#blurb = "some text"
#ax1.text(0.5, 0.5, blurb, transform=ax1.transAxes)
#ax1.text(0.5, 0.5, blurb, transform=ax1.transAxes,
#      va='top', ha='left', bbox=dict(facecolor='white', pad=10.0))
#      fontdict={'family':'monospace'}) # fixed-width

fig.tight_layout() # adjust boundaries sensibly, matplotlib v1.1+
plt.draw()
if context.plot_name:
    sys.stdout.write("Saving diagram to %s ... " % context.plot_name)
    fig.savefig(context.plot_name, bbox_inches='tight')
    sys.stdout.write("done.\n")

# cyclical colormap ... cmocean.cm.phase
# cmocean: https://matplotlib.org/cmocean/




######################################################################
# CHANGELOG (kwhere.py):
#---------------------------------------------------------------------
#
#  2018-06-01:
#     -- Increased __version__ to 0.2.1.
#     -- Switched RA/Dec labels to +/- formatting.
#
#  2018-05-23:
#     -- Increased __version__ to 0.2.0.
#     -- Useful plotting functionality now provided. Script cleaned up.
#
#  2018-05-22:
#     -- Increased __version__ to 0.1.0.
#     -- First created kwhere.py.
#
