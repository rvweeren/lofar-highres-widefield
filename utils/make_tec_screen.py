from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.pyplot import figure
from scipy.interpolate import Rbf

import casacore.tables as ct
import concurrent.futures
import losoto.h5parm as h5parm
import matplotlib.pyplot as plt
import numpy as np

import sys


def plot_screen(img, prefix='', title='', suffix='', wcs=None):
    fig = figure(figsize=(12, 12))
    if wcs is not None:
        ax = fig.add_subplot(111, projection=wcs)
    else:
        ax = fig.add_subplot(111, projection=wcs)

    ax.set_title(title)
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')
    i = ax.imshow(img, origin='lower', vmin=-0.1, vmax=0.1)
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
if np.rad2deg(phasecenter[0]) < 0.: # avoid negative CRVAL1
   phasecenter[0] = phasecenter[0] + (2.*np.pi)

t.close()
h5 = h5parm.h5parm(h5p)
ss = h5.getSolset('sol000')
st = ss.getSoltab('tec000')
time = st.getAxisValues('time')
stime = time[0]
dtime = time[1] - time[0]

Nantenna = len(names)
Ntimes = len(time)  # 60s timeslots for an 8 hour pointing
# Set the frequency axis (that TEC doens't use) to 150 MHz as central frequency with 50 MHz of bandwidth.
header='''SIMPLE  =                    T / file does conform to FITS standard
BITPIX  =                  -32 / number of bits per data pixel
NAXIS   =                    5 / number of data axes
NAXIS1  =                 256 / length of RA axis
NAXIS2  =                 256 / length of DEC axis
NAXIS3  =                   {nant:d} / length of ANTENNA axis
NAXIS4  =                    1 / length of FREQ axis
NAXIS5  =                   {ntimes:d} / length of TIME axis
EXTEND  =                    T / FITS dataset may contain extensions
EQUINOX =                2000.
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
CDELT5  =                  {dtime:f}'''.format(nant=Nantenna, ntimes=Ntimes, ra=np.rad2deg(phasecenter[0]), dec=np.rad2deg(phasecenter[1]), cfreq=150e6, dfreq=100e6, stime=stime, dtime=dtime)

# 4ch/4s data
# Frequency has length one, since it is TEC.
# FITS is transposed compared to Numpy.
data = np.zeros((Ntimes, 1, Nantenna, 256, 256), dtype=np.float32)

# Read in h5parm.
h5 = h5parm.h5parm(h5p)
ss = h5.getSolset('sol000')
st = ss.getSoltab('tec000')
h5_stations = list(st.getAxisValues('ant'))
# Find nearest pixel for a given direction.
H = fits.Header.fromstring(header, sep='\n')
wcs = WCS(H).celestial
directions = list(st.getAxisValues('dir'))
print(directions)
sources = ss.getSou()
RA = []
DEC = []
TEC = st.getValues()[0]
dirs = []
print(TEC.shape)
for d in directions:
    c = sources[d]
    diridx = directions.index(d)
    dirs.append(d)
    RAd, DECd = np.rad2deg(c)
    RA.append(RAd)
    DEC.append(DECd)
    c = wcs.wcs_world2pix(RAd, DECd, 0)
    print(c[0])
    print('Direction {:s} is at pixel coordinates {:f},{:f}'.format(d, float(c[0]), float(c[1])))

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
            idx = h5_stations.index('ST001')
        else:
            idx = h5_stations.index(station)
        refidx = h5_stations.index('ST001')
        vals_tmp = st.getValues()[0][diridx, idx, :] - st.getValues()[0][diridx, refidx, :]
        data[:, :, istation, Y, X] = vals_tmp.reshape((-1, 1))

hdu = fits.PrimaryHDU(header=H)
hdu.data = data
hdu.writeto('tecscreen_raw.fits')
RA = np.asarray(RA)
DEC = np.asarray(DEC)
# Interpolate the grid using a nearest neighbour approach.
# https://stackoverflow.com/questions/5551286/filling-gaps-in-a-numpy-array


def interpolate_station(antname):
    print('Processing antenna {:s}.'.format(antname))
    refidx = h5_stations.index('ST001')
    if 'CS' in antname:
        # Take ST001 solutions.
        interpidx = h5_stations.index('ST001')
    else:
        interpidx = h5_stations.index(antname)
    # These will be the new coordinates to evaluate the screen over.
    size = 256
    # Do the interpolation in pixel space.
    x = np.arange(size)
    y = np.arange(size)
    xx, yy = np.meshgrid(x, y)
    # Do radial basis function interpolation for each time step.
    screen = np.zeros((TEC.shape[-1], 256, 256))
    # Iterate over all timeslots.
    for itstep in range(TEC.shape[-1]):
        X, Y = np.around(wcs.wcs_world2pix(RA, DEC, 0)).astype(int)
        tecs = TEC[:, interpidx, itstep] - TEC[:, refidx, itstep]
        rbf = Rbf(X, Y, tecs, smooth=0.0)
        '''
        if ('CS' not in antname) and ('RS' not in antname):
            print('== bla ==')
            for i, (rr,dd) in enumerate(zip(RA,DEC)):
                print(rr)
                print(dd)
                orig = tstep[interpidx][i]
                interp = rbf(rr, dd)
                print('Difference TEC and Rbf: {:e}'.format(orig - interp))
                print('== end bla ==')
        '''
        tinterp = rbf(xx, yy)
        if not np.allclose(tinterp[Y, X], tecs, rtol=1e-5):
            raise ValueError('Interpolated screen does not go through nodal points.')
        screen[itstep, :, :] = tinterp
    '''
    didx = directions.index('P470')
    print(didx)
    print(tecs.shape)
    print(X[didx], Y[didx])
    print(screen[:, 136, 129].shape)
    fig = figure()
    fig.suptitle(antname)
    ax = fig.add_subplot(111)
    ax.plot(TEC[didx, interpidx, :] - TEC[didx, refidx, :], 'C0h--', label='H5parm', markersize=12)
    ax.plot(screen[:, 136, 129], 'C1h--', label='Rbf Interp. Y,X')
    ax.plot(screen[:, 129, 136], 'C2h--', label='Rbf Interp. X,Y')
    ax.legend()
    ax.set_xlabel('Time'); ax.set_ylabel('dTEC')
    show()
    del fig
    if ('CS' not in antname) and ('RS' not in antname):
        print('== bla ==')
        for i, (rr,dd) in enumerate(zip(RA,DEC)):
            print(rr)
            print(dd)
            orig = tstep[interpidx][i]
            interp = rbf(rr, dd)
            print('Difference TEC and Rbf: {:e}'.format(orig - interp))
            print('== end bla ==')
    fig = figure()
    fig.suptitle(antname)
    ax = fig.add_subplot(111)
    im = ax.imshow(screen[..., -1], origin='lower', extent=[ra_max, ra_min, dec_min, dec_max])
    ax.scatter(RA, DEC, marker='x', color='r')
    ax.set_xlabel('Right ascension'); ax.set_ylabel('Declination')
    fig.colorbar(im)
    fig.savefig(antname + '_{:03d}.png'.format(itstep))
    del fig
    show()
    '''
    return names.index(antname), screen


# Inspired by https://stackoverflow.com/questions/38309535/populate-numpy-array-through-concurrent-futures-multiprocessing
data_int = np.zeros(data.shape)
print('Making TEC screen.')
# Run each antenna in parallel.
with concurrent.futures.ProcessPoolExecutor(64) as executor:
    for antenna, screen in executor.map(interpolate_station, names):
        screen = screen
        data_int[:, 0, antenna, :, :] = screen

hdu = fits.PrimaryHDU(header=H)
hdu.data = data_int
hdu.writeto('tecscreen_rbf.fits')
print('Finished interpolating.')
