from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage

import casacore.tables as ct
import losoto.h5parm as h5parm
import matplotlib.pyplot as plt
import numpy as np

import sys


def plot_screen(img, prefix='', title='', suffix='', wcs=None):
    fig = figure(figsize=(12,12))
    if wcs is not None:
        ax = fig.add_subplot(111, projection=wcs)
    else:
        ax = fig.add_subplot(111, projection=wcs)
    
    ax.set_title(title)
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')
    i = ax.imshow(img, origin='lower', vmin=-0.075, vmax=0.075)
    cbar = fig.colorbar(i, fraction=0.046, pad=0.04)
    cbar.set_label('dTEC')
    fig.savefig(prefix + 'tecscreen_' + suffix + '.svg')
    plt.close(fig)
    del fig


ms = sys.argv[1]
h5p = sys.argv[2]

t = ct.taql('SELECT NAME FROM {ms:s}::ANTENNA'.format(ms=ms))
names = t.getcol('NAME')
t.close()
t = ct.taql('SELECT REFERENCE_DIR FROM {ms:s}::FIELD'.format(ms=ms))
phasecenter = t.getcol('REFERENCE_DIR').squeeze()
print(phasecenter)
t.close()
# Time is stored as MJD, convert from seconds to days here as that's what FITS wants.
t = ct.taql('SELECT TIME FROM {ms:s}'.format(ms=ms))
time = t.getcol('TIME')
stime = time / (24 * 60 * 60)
dtime = (time[1] - time[0]) / (24 * 60 * 60)
t.close()
Nantenna = len(names)
Ntimes = 480 # 60s timeslots for an 8 hour pointing
# Pull values from the measurement set normally.
# Pull values from the measurement set normally.
header='''SIMPLE  =                    T / file does conform to FITS standard
BITPIX  =                  -32 / number of bits per data pixel
NAXIS   =                    5 / number of data axes
NAXIS1  =                 128 / length of RA axis
NAXIS2  =                 128 / length of DEC axis
NAXIS3  =                   {nant:d} / length of ANTENNA axis
NAXIS4  =                    1 / length of FREQ axis
NAXIS5  =                   {ntimes:d} / length of TIME axis
EXTEND  =                    T / FITS dataset may contain extensions
CTYPE1  = 'RA---SIN'           / Right ascension angle cosine
CRPIX1  =                 65.
CRVAL1  =          {ra:f}
CDELT1  =              -0.0195
CUNIT1  = 'deg     '
CTYPE2  = 'DEC--SIN'           / Declination angle cosine
CRPIX2  =                 65.
CRVAL2  =     {dec:f}
CDELT2  =               0.0195
CUNIT2  = 'deg     '
CTYPE3  = 'ANTENNA '
CRPIX3  =                   1.
CRVAL3  =                   0.
CTYPE4  = 'FREQ    '           / Central frequency
CRPIX4  =                   1.
CRVAL4  =     {cfreq:f}
CDELT4  =         {dfreq:f}
CUNIT4  = 'Hz      '
CTYPE5  = 'TIME    '
CRPIX5  =                   1.
CRVAL5  =                   {stime:f} / Should be an AIPS time
CDELT5  =                  {dtime:f}'''.format(nant=Nantenna, ntimes=Ntimes, ra=np.rad2deg(phasecenter[0]), dec=np.rad2deg(phasecenter[1]), cfreq=144627380.37109, dfreq=48828.1, stime=5038110492.002781, dtime=60.)

# 4ch/4s data
# Frequency has length one, since it is TEC.
#data = np.zeros((128, 128, Nantenna, 1, Ntimes))
# FITS is transposed compared to Numpy.
data = np.zeros((Ntimes, 1, Nantenna, 128, 128))

# Read in h5parm.
h5 = h5parm.h5parm(h5p)
ss = h5.getSolset('sol000')
st = ss.getSoltab('tec000')
h5_stations = list(st.getAxisValues('ant'))
# Find nearest pixel for a given direction.
H = fits.Header.fromstring(header, sep='\n')
wcs = WCS(H).celestial
directions = list(st.getAxisValues('dir'))
for d, c in ss.getSou().items():
    diridx = directions.index(d)
    RAd, DECd = np.rad2deg(c)
    print('Adding direction {:s} at {:f},{:f}'.format(d, RAd, DECd))
    pixel_coords = wcs.wcs_world2pix(RAd, DECd, 0)
    X = int(pixel_coords[0])
    Y = int(pixel_coords[1])
    print('== Direction corresponds to (x, y) = ({:d}, {:d})'.format(X, Y))
    # Fill it with the solutions for each antenna.
    for istation, station in enumerate(names):
        if 'CS' in station:
            # Take ST001 solutions.
            # [dir, ant, time]
            idx = h5_stations.index('ST001')
        else:
            idx = h5_stations.index(station)
        vals_tmp = st.getValues()[0][diridx, idx, :].reshape((-1, 1))
        data[:, :, istation, Y, X] = vals_tmp

hdu = fits.PrimaryHDU(header=H)
hdu.data = data
hdu.writeto('tecscreen_raw.fits')

def interpolate_antenna(data, antenna):
    pass

# Interpolate the grid using a nearest neighbour approach.
# https://stackoverflow.com/questions/5551286/filling-gaps-in-a-numpy-array
import SharedArray as sa
#data_int = sa.create('shm://interpdata', data.shape)
#print('Interpolating from given directions using nearest neighbour.')

import multiprocessing as mp
output = mp.Queue()
from matplotlib.pyplot import figure

data_int = np.zeros(data.shape)
print('Making TEC screen.')
# Run the antennas in parallel. Much faster than interpolating an NxNxNantx1xNtime array.
for i, a in enumerate(names):
    print('Processing antenna {:s}'.format(a))
    adata = data[:, :, i, :, :]
    invalid_cell_mask = np.where(adata == 0, 1, 0)
    indices = ndimage.distance_transform_edt(invalid_cell_mask, return_distances=False, return_indices=True)
    data_int[:, :, i, :, :] = adata[tuple(indices)]
    #for j in range(Ntimes):
    #    plot_screen(data_int[j, 0, i, :, :], title=a, prefix=a, suffix='{:03d}'.format(j), wcs=wcs)
hdu = fits.PrimaryHDU(header=H)
hdu.data = data_int
hdu.writeto('tecscreen_interpolated.fits')
print(data_int.shape)
print('Finished interpolating.')

