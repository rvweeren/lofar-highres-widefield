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
parser.add_argument('--convert-tec', dest='convert_tec', action='store_true', default=False)
args = parser.parse_args()
convert_tec = args.convert_tec
print(convert_tec)

mslist = sorted(glob.glob(args.msdir + '/*.' + args.mssuffix))
ms_first = mslist[0]
ms_last = mslist[-1]
print(ms_first, ms_last)

print('Determining time grid...')
'''
tf = ct.taql('SELECT UNIQUE TIME FROM $ms_first')
time = tf.getcol('TIME')
time_first = time[0]
time_last = time[-1]
time_spacing = time[1] - time[0]
time_spacing = 60.
tf.close()
print(time_first, time_last, time_spacing)
'''


#ax_time = np.arange(time_first, time_last + time_spacing, time_spacing)
h5 = h5parm.h5parm('dummy.h5')
ss = h5.getSolset('sol000')
st = ss.getSoltab('tec000')
ax_time = st.getAxisValues('time')

#time = st.getAxisValues('time')
#freq = st.getAxisValues('freq')
#time_first = time[0]; time_last = time[-1]; time_spacing = time[1] - time[0]

h5 = h5parm.h5parm(args.h5parms[0])
ss = h5.getSolset(args.solsetin)
st = ss.getSoltab(args.soltab2merge)

#freq_first = freq[0]; freq_last = freq[-1]; freq_spacing = freq[1] - freq[0]

#ax_time = np.arange(time_first, time_last + time_spacing, time_spacing)

antennas = st.getAxisValues('ant')
ss_antennas = ss.obj.antenna.read()
directions = []

vals = st.getValues()[0]
#vals_reordered = reorderAxes(vals, st.getAxesNames(), ['time', 'freq', 'ant', 'dir'])
vals_reordered = reorderAxes(vals, st.getAxesNames(), ['dir', 'ant', 'freq', 'time'])
#phases = np.zeros((len(ax_time), len(ax_freq), len(antennas), len(args.h5parms)))
#phases = np.zeros((1, len(args.h5parms), len(antennas), len(ax_freq), len(ax_time)))
if convert_tec:
    print('Determining frequency grid...')
    ff = ct.taql('SELECT CHAN_FREQ, CHAN_WIDTH FROM ' + ms_first + '::SPECTRAL_WINDOW')
    freq_first = ff.getcol('CHAN_FREQ')[0][0]
    freq_spacing = ff.getcol('CHAN_WIDTH')[0][0]
    ff.close()

    fl = ct.taql('SELECT CHAN_FREQ, CHAN_WIDTH FROM ' + ms_last + '::SPECTRAL_WINDOW')
    freq_last = fl.getcol('CHAN_FREQ')[0][0]
    print(freq_first, freq_last, freq_spacing)
    ax_freq = np.arange(freq_first, freq_last + freq_spacing, freq_spacing)
    phases = np.zeros((1, 1, len(antennas), len(ax_freq), len(ax_time)))
elif not convert_tec:
    phases = np.zeros((1, len(antennas), len(ax_time)))
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
    source_coords = ss.getSou()['pointing']
    d = 'Dir{:02d}'.format(i)
    if (d, source_coords) not in sourcelist:
        print('Adding new direction {:f},{:f}'.format(*source_coords))
        idx = i
        directions.append(d)
        sourcelist.append((d, source_coords))
    else:
        # Direction already exists, add to the existing solutions.
        print('Direction {:f},{:f} already exists, adding solutions instead.'.format(*source_coords))
        idx = directions.index(d)
        d = 'Dir{:02d}'.format(idx)
    print(st.getAxisValues('dir'))
    if st.getType() == 'tec' and convert_tec:
        # Convert tec to phase.
        tec_tmp = st.getValues()[0]
        #tec = reorderAxes(tec_tmp, st.getAxesNames(), ['time', 'freq', 'ant', 'dir'])
        tec = reorderAxes(tec_tmp, st.getAxesNames(), ['dir', 'ant', 'freq', 'time'])
        # -1 assumes the expected shape along the frequency axis.
        #freqs = ax_freq.reshape(1, -1, 1, 1, 1)
        #freqs = ax_freq.reshape(1, 1, -1, 1)
        freqs = ax_freq.reshape(1, 1, 1, -1, 1)
        tecphase = (-8.4479745e9 * tec / freqs)
        tp = interp_along_axis(tecphase, st.getAxisValues('time'), ax_time, -1)
        # Now add the phases to the total phase correction for this direction.
        #phases[:, :, :, i] += tp[..., 0]
        #phases[0, idx, :, :, :] += tp[0, ...]
        if idx == 0:
            phases[:, idx, :, :, :] += tp[:, 0, :, :, :]
        elif idx > 0:
            phases = np.append(phases, tp, axis=1)
    elif st.getType() == 'tec' and not convert_tec:
        tec_tmp = st.getValues()[0]
        tec = reorderAxes(tec_tmp, st.getAxesNames(), ['dir', 'ant', 'freq', 'time'])[0, :, 0, :]
        # -1 assumes the expected shape along the frequency axis.
        tp = interp_along_axis(tec, st.getAxisValues('time'), ax_time, -1)
        tp = tp.reshape(-1, *tp.shape)
        print(tp.shape)
        print(phases.shape)
        # Now add the tecs to the total phase correction for this direction.
        if idx == 0:
            # Axis order is pol,dir,ant,time.
            # Set the first direction.
            print(tp.shape)
            print(phases[idx, :, :].shape)
            phases[idx, :, :] += tp[0, :, :]
        elif idx > 0:
            phases = np.append(phases, tp, axis=0)
    elif st.getType() == 'phase':
        # Convert tec to phase.
        phase_tmp = st.getValues()[0]
        #tec = reorderAxes(tec_tmp, st.getAxesNames(), ['time', 'freq', 'ant', 'dir'])
        phase = reorderAxes(phase_tmp, st.getAxesNames(), ['dir', 'ant', 'freq', 'time'])
        # -1 assumes the expected shape along the frequency axis.
        #freqs = ax_freq.reshape(1, -1, 1, 1, 1)
        tp = interp_along_axis(phase, st.getAxisValues('time'), ax_time, -1)
        # Now add the phases to the total phase correction for this direction.
        #phases[:, :, :, i] += tp[..., 0]
        phases[0, idx, :, :, :] += tp[0, ...]
    h5.close()

# Create the output h5parm.
weights = np.ones(phases.shape)
sources = np.array(sourcelist, dtype=[('name', 'S128'), ('dir', '<f4', (2,))])
solsetout.obj.source.append(sources)
if not convert_tec:
    solsetout.makeSoltab('tec', axesNames=['dir', 'ant', 'time'], axesVals=[directions, antennas, ax_time], vals=phases, weights=weights)
elif convert_tec:
    solsetout.makeSoltab('phase', axesNames=['pol', 'dir', 'ant', 'freq', 'time'], axesVals=[['XX'], directions, antennas, ax_freq, ax_time], vals=phases, weights=weights)
h5out.close()
