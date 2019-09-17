from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from losoto.lib_operations import reorderAxes
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
    i = ax.imshow(img, origin='lower', vmin=0.5, vmax=1.5)
    cbar = fig.colorbar(i, fraction=0.046, pad=0.04)
    cbar.set_label('Gain')
    fig.savefig(prefix + 'gainscreen_' + suffix + '.svg')
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

h5 = h5parm.h5parm(h5p)
ss = h5.getSolset('sol000')
st = ss.getSoltab('amplitude000')
time = st.getAxisValues('time')
stime = time[0] 
dtime = time[1] - time[0]
freq = st.getAxisValues('freq')
sfreq = freq[0]
dfreq = freq[1] - freq[0]

Nantenna = len(names)
Ntimes = len(time) # 60s timeslots for an 8 hour pointing
Nfreqs = len(freq)
# Set the frequency axis (that TEC doens't use) to 150 MHz as central frequency with 50 MHz of bandwidth.
header='''SIMPLE  =                    T / file does conform to FITS standard
BITPIX  =                  -32 / number of bits per data pixel
NAXIS   =                    6 / number of data axes
NAXIS1  =                 256 / length of RA axis
NAXIS2  =                 256 / length of DEC axis
NAXIS3  =                   4
NAXIS4  =                   {nant:d} / length of ANTENNA axis
NAXIS5  =                    1 / length of FREQ axis
NAXIS6  =                   {ntimes:d} / length of TIME axis
EXTEND  =                    T / FITS dataset may contain extensions
CTYPE1  = 'RA---SIN'           / Right ascension angle cosine
CRPIX1  =                 128.
CRVAL1  =          {ra:f}
CDELT1  =              -0.0195
CUNIT1  = 'deg     '
CTYPE2  = 'DEC--SIN'           / Declination angle cosine
CRPIX2  =                 128.
CRVAL2  =     {dec:f}
CDELT2  =               0.0195
CUNIT2  = 'deg     '
CTYPE3  = 'MATRIX  '
CRPIX3  =                   1.
CDELT3  =                   1.
CTYPE4  = 'ANTENNA '
CRPIX4  =                   1.
CRVAL4  =                   0.
CTYPE5  = 'FREQ    '           / Central frequency
CRPIX5  =                   1.
CRVAL5  =     {cfreq:f}
CDELT5  =         {dfreq:f}
CUNIT5  = 'Hz      '
CTYPE6  = 'TIME    '
CRPIX6  =                   1.
CRVAL6  =                   {stime:f} / Should be an AIPS time
CDELT6  =                  {dtime:f}'''.format(nant=Nantenna, ntimes=Ntimes, ra=np.rad2deg(phasecenter[0]), dec=np.rad2deg(phasecenter[1]), cfreq=sfreq, dfreq=dfreq, stime=stime, dtime=dtime)

# 4ch/4s data
# FITS is transposed compared to Numpy.
# 4 is the "matrix" entry containing Re(XX), Im(XX), Re(YY), Im(YY)
data = np.zeros((Ntimes, Nfreqs, Nantenna, 4, 256, 256), dtype=np.float32)


# Read in h5parm.
h5 = h5parm.h5parm(h5p)
ss = h5.getSolset('sol000')
sta = ss.getSoltab('amplitude000')
stp = ss.getSoltab('phase000')
h5_stations = list(st.getAxisValues('ant'))
H = fits.Header.fromstring(header, sep='\n')
wcs = WCS(H).celestial
directions = list(st.getAxisValues('dir'))
# Find nearest pixel for a given direction.
for d, c in ss.getSou().items():
    diridx = directions.index(d)
    RAd, DECd = np.rad2deg(c)
    print('Adding direction {:s} at {:f},{:f}'.format(d, RAd, DECd))
    pixel_coords = wcs.wcs_world2pix(RAd, DECd, 0)
    X = int(np.around(pixel_coords[0]))
    Y = int(np.around(pixel_coords[1]))
    print('== Direction corresponds to (x, y) = ({:d}, {:d})'.format(X, Y))
    # Fill it with the solutions for each antenna.
    for istation, station in enumerate(names):
        if 'CS' in station:
            # Take ST001 solutions.
            # [dir, ant, time]
            sidx = h5_stations.index('ST001')
        else:
            sidx = h5_stations.index(station)
        valsa_tmp = sta.getValues()[0]
        valsa_tmp_r = reorderAxes(valsa_tmp, sta.getAxesNames(), ['time', 'freq', 'ant', 'dir', 'pol'])

        valsp_tmp = stp.getValues()[0]
        valsp_tmp_r = reorderAxes(valsp_tmp, stp.getAxesNames(), ['time', 'freq', 'ant', 'dir', 'pol'])

        real = valsa_tmp_r * np.cos(valsp_tmp_r)
        imag = valsa_tmp_r * np.sin(valsp_tmp_r)

        reXX = real[..., 0]
        imXX = imag[..., 0]
        reYY = real[..., -1]
        imYY = imag[..., -1]
        matrix = np.asarray([reXX, imXX, reYY, imYY], dtype=np.float32)
        #print('Matrix: ', matrix.shape)
        matrix = np.transpose(matrix, (1, 2, 3, 0, 4))
        #print('Matrix.T: ', matrix.shape)
        #print('Data: ', data.shape)
        data[:, :, istation, :, Y, X] = matrix[:, :, sidx, :, diridx]

hdu = fits.PrimaryHDU(header=H)
hdu.data = data
hdu.writeto('gainscreen_raw_40dirs_bigger.fits', overwrite=True)

# Interpolate the grid using a nearest neighbour approach.
# https://stackoverflow.com/questions/5551286/filling-gaps-in-a-numpy-array
import SharedArray as sa
#data_int = sa.create('shm://interpdata', data.shape)
#print('Interpolating from given directions using nearest neighbour.')

import multiprocessing as mp
output = mp.Queue()
from matplotlib.pyplot import figure

data_int = np.zeros(data.shape)
print('Making gain screen.')
# Run the antennas in parallel. Much faster than interpolating an NxNxNantx1xNtime array.
for i, a in enumerate(names):
    print('Processing antenna {:s}'.format(a))
    adata = data[:, :, i, :, :, :]
    invalid_cell_mask = np.where(adata == 0, 1, 0)
    indices = ndimage.distance_transform_edt(invalid_cell_mask, return_distances=False, return_indices=True)
    data_int[:, :, i, :, :, :] = adata[tuple(indices)]
    if 'CS' not in a:
        for j in range(Ntimes):
            real = data_int[j, 0, i, 0, :, :]
            imag = data_int[j, 0, i, 1, :, :]
            amp = (real**2 + imag**2)**0.5
            #plot_screen(amp, title=a + ' XX', prefix=a, suffix='{:03d}'.format(j), wcs=wcs)
# We don't want 0 gains, so check for zeros.
data_int[np.where(data_int == 0)] = 1
hdu = fits.PrimaryHDU(header=H)
hdu.data = data_int
hdu.writeto('gainscreen_interpolated_40dirs_bigger.fits', overwrite=True)
print('Finished interpolating.')

