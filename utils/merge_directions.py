import argparse
import glob
import sys

from losoto.lib_operations import reorderAxes
from scipy.interpolate import interp1d

import casacore.tables as ct
import losoto.h5parm as h5parm
import numpy as np


def interp_along_axis(x, interp_from, interp_to, axis):
    print('Inter/extrapolating from {:d} to {:d} along axis {:d}'.format(len(interp_from), len(interp_to), axis))
    interp_vals = interp1d(interp_from, x, axis=axis, kind='nearest', fill_value='extrapolate')
    new_vals = interp_vals(interp_to)
    return new_vals

parser =argparse.ArgumentParser()
parser.add_argument('--mspath', dest='msdir')
parser.add_argument('--mssuffix', dest='mssuffix', default='ms')
parser.add_argument('--h5parms', dest='h5parms', nargs='+')
parser.add_argument('--soltab', dest='soltab2merge')
parser.add_argument('--solset-in', dest='solsetin')
parser.add_argument('--h5parm-out', dest='h5out')
args = parser.parse_args()

mslist = sorted(glob.glob(args.msdir + '/*.' + args.mssuffix))
ms_first = mslist[0]
ms_last = mslist[-1]
print(ms_first, ms_last)

print('Determining time grid...')
tf = ct.taql('SELECT UNIQUE TIME FROM $ms_first')
time = tf.getcol('TIME')
time_first = time[0]
time_last = time[-1]
time_spacing = time[1] - time[0]
tf.close()
print(time_first, time_last, time_spacing)

print('Determining frequency grid...')
ff = ct.taql('SELECT CHAN_FREQ, CHAN_WIDTH FROM ' + ms_first + '::SPECTRAL_WINDOW')
freq_first = ff.getcol('CHAN_FREQ')[0][0]
freq_spacing = ff.getcol('CHAN_WIDTH')[0][0]
ff.close()

fl = ct.taql('SELECT CHAN_FREQ, CHAN_WIDTH FROM ' + ms_last + '::SPECTRAL_WINDOW')
freq_last = fl.getcol('CHAN_FREQ')[0][0]
print(freq_first, freq_last, freq_spacing)

ax_time = np.arange(time_first, time_last + time_spacing, time_spacing)
ax_freq = np.arange(freq_first, freq_last + freq_spacing, freq_spacing)

h5 = h5parm.h5parm(args.h5parms[0])
ss = h5.getSolset(args.solsetin)
st = ss.getSoltab(args.soltab2merge)
antennas = st.getAxisValues('ant')
ss_antennas = ss.obj.antenna.read()
directions = []

vals = st.getValues()[0]
#vals_reordered = reorderAxes(vals, st.getAxesNames(), ['time', 'freq', 'ant', 'dir'])
vals_reordered = reorderAxes(vals, st.getAxesNames(), ['dir', 'ant', 'freq', 'time'])
#phases = np.zeros((len(ax_time), len(ax_freq), len(antennas), len(args.h5parms)))
phases = np.zeros((1, len(args.h5parms), len(antennas), len(ax_freq), len(ax_time)))
#print('Vals shape: ', vals.shape)
#print('Reordered shape: ', vals_reordered.shape)
#print('Values shape: ', values.shape)
h5out = h5parm.h5parm(args.h5out, readonly=False)
# Create a new source table, ripped from the LoSoTo code.
#descriptor = np.dtype([('name', np.str_, 128),('dir', np.float32, 2)])
#sources = h5out.H.create_table(solset, 'source', descriptor, title = 'Source names and directions', expectedrows = 25)
solsetout = h5out.makeSolset('sol000')
antennasout = solsetout.getAnt()
#solsetout.obj.antenna.append(ss.obj.antenna.col('name'))
antennatable = solsetout.obj._f_get_child('antenna')
antennatable.append(ss.obj.antenna.read())
sourcelist = []

for i, h5 in enumerate(args.h5parms):
    print('Processing direction for ' + h5)
    # Read in the data
    h5 = h5parm.h5parm(h5)
    ss = h5.getSolset(args.solsetin)
    st = ss.getSoltab(args.soltab2merge)
    directions.append('Dir{:02d}'.format(i))
    print(st.getAxisValues('dir'))
    if st.getType() == 'tec':
        # Convert tec to phase.
        tec_tmp = st.getValues()[0]
        #tec = reorderAxes(tec_tmp, st.getAxesNames(), ['time', 'freq', 'ant', 'dir'])
        tec = reorderAxes(tec_tmp, st.getAxesNames(), ['dir', 'ant', 'freq', 'time'])
        # -1 assumes the expected shape along the frequency axis.
        #freqs = ax_freq.reshape(1, -1, 1, 1, 1)
        freqs = ax_freq.reshape(1, 1, -1, 1)
        tecphase = (-8.4479745e9 * tec / freqs)
        tp = interp_along_axis(tecphase, st.getAxisValues('time'), ax_time, -1)
        # Now add the phases to the total phase correction for this direction.
        #phases[:, :, :, i] += tp[..., 0]
        phases[0, i, :, :, :] += tp[0, ...]
    source_coords = ss.getSou()['pointing']
    sourcelist.append(('Dir{:02d}'.format(i), source_coords))
    h5.close()

# Create the output h5parm.
weights = np.ones(phases.shape)
sources = np.array(sourcelist, dtype=[('name', 'S128'), ('dir', '<f4', (2,))])
solsetout.obj.source.append(sources)
#solsetout.makeSoltab('phase', axesNames=['time', 'freq', 'ant', 'dir'], axesVals=[ax_time, ax_freq, antennas, directions], vals=phases, weights=weights)
solsetout.makeSoltab('phase', axesNames=['pol', 'dir', 'ant', 'freq', 'time'], axesVals=[['XX'], directions, antennas, ax_freq, ax_time], vals=phases, weights=weights)
h5out.close()
