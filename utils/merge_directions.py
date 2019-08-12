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
parser.add_argument('--mspath', dest='msdir', description='Path to the directory with easurement sets to pull frequency axis from, when converting TEC to phase.')
parser.add_argument('--mssuffix', dest='mssuffix', default='ms', description='Suffix of your measurement sets, e.g. MS or ms.')
parser.add_argument('--h5parms', dest='h5parms', nargs='+', description='Input H5parms to merge as directions, where each h5parm is one direction.')
parser.add_argument('--soltab', dest='soltab2merge', description='SolTab of the H5parms to merge.')
parser.add_argument('--solset-in', dest='solsetin', description='SolSet to take the soltab from.')
parser.add_argument('--h5parm-out', dest='h5out', description='Output H5parm with all directions present.')
parser.add_argument('--convert-tec', dest='convert_tec', action='store_true', default=False, description='Convert TEC values to their corresponding phase corrections base on the frequencies in the Measurement Sets.')
args = parser.parse_args()
convert_tec = args.convert_tec

mslist = sorted(glob.glob(args.msdir + '/*.' + args.mssuffix))
ms_first = mslist[0]
ms_last = mslist[-1]

h5 = h5parm.h5parm('dummy.h5')
ss = h5.getSolset('sol000')
st = ss.getSoltab('tec000')
print('Determining time grid...')
ax_time = st.getAxisValues('time')
h5.close()

h5 = h5parm.h5parm(args.h5parms[0])
ss = h5.getSolset(args.solsetin)
st = ss.getSoltab(args.soltab2merge)

antennas = st.getAxisValues('ant')
ss_antennas = ss.obj.antenna.read()
directions = []

vals = st.getValues()[0]
AN = st.getAxesNames()
axes_new = ['ant', 'freq', 'time']
if 'dir' in AN:
    axes_new = ['dir'] + axes_new
if 'pol' in AN:
    axes_new = ['pol'] + axes_new
    polarizations = st.getAxisValues('pol')

vals_reordered = reorderAxes(vals, st.getAxesNames(), axes_new)

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
    phases = np.zeros(vals_reordered.shape)

h5out = h5parm.h5parm(args.h5out, readonly=False)
solsetout = h5out.makeSolset('sol000')
antennasout = solsetout.getAnt()
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
    if st.getType() == 'tec' and convert_tec:
        # Convert tec to phase.
        tec_tmp = st.getValues()[0]
        #tec = reorderAxes(tec_tmp, st.getAxesNames(), ['time', 'freq', 'ant', 'dir'])
        tec = reorderAxes(tec_tmp, st.getAxesNames(), axes_new)
        # -1 assumes the expected shape along the frequency axis.
        freqs = ax_freq.reshape(1, 1, 1, -1, 1)
        tecphase = (-8.4479745e9 * tec / freqs)
        tp = interp_along_axis(tecphase, st.getAxisValues('time'), ax_time, -1)
        # Now add the phases to the total phase correction for this direction.
        if idx == 0:
            phases[:, idx, :, :, :] += tp[:, 0, :, :, :]
        elif idx > 0:
            phases = np.append(phases, tp, axis=1)
    elif st.getType() == 'tec' and not convert_tec:
        tec_tmp = st.getValues()[0]
        if 'dir' in axes_new:
            # TEC will never have a polarization axis and has only one frequency, so the latter indexing can stay: first direction, all antennas, first frequency, all times.
            tec = reorderAxes(tec_tmp, st.getAxesNames(), axes_new)[0, :, 0, :]
        else:
            tec = reorderAxes(tec_tmp, st.getAxesNames(), axes_new)
        # -1 assumes the expected shape along the frequency axis.
        tp = interp_along_axis(tec, st.getAxisValues('time'), ax_time, -1)
        tp = tp.reshape(-1, *tp.shape)
        # Now add the tecs to the total phase correction for this direction.
        if idx == 0:
            # Axis order is dir,ant,time.
            # Set the first direction.
            if 'dir' in axes_new:
                phases[idx, :, :] += tp[0, :, :]
            else:
                phases[idx, :, :] = tp
        elif idx > 0:
            phases = np.append(phases, tp, axis=0)
    elif st.getType() == 'phase':
        phase_tmp = st.getValues()[0]
        phase = reorderAxes(phase_tmp, st.getAxesNames(), axes_new)
        tp = interp_along_axis(phase, st.getAxisValues('time'), ax_time, -1)
        # Now add the phases to the total phase correction for this direction.
        if pol in axes_new:
            polidx = range(phases.shape[0])
        else:
            polidx = 0
        for pidx in polidx:
            if 'dir' in axes_new:
                phases[pidx, idx, :, :, :] += tp[0, ...]
            else:
                phases[pidx, idx, :, :, :] += tp
    h5.close()

# Create the output h5parm.
weights = np.ones(phases.shape)
sources = np.array(sourcelist, dtype=[('name', 'S128'), ('dir', '<f4', (2,))])
solsetout.obj.source.append(sources)
if args.soltab2merge == 'phase' and len(polarizations) > 0:
    solsetout.makeSoltab('phase', axesNames=axes_new, axesVals=[polarizations, directions, antennas, ax_freq, ax_time], vals=phases, weights=weights)
elif args.soltab2merge == 'phase' and len(polarizations) == 0:
    solsetout.makeSoltab('phase', axesNames=axes_new, axesVals=[directions, antennas, ax_freq, ax_time], vals=phases[0,...], weights=weights)
if not convert_tec:
    solsetout.makeSoltab('tec', axesNames=['dir', 'ant', 'time'], axesVals=[directions, antennas, ax_time], vals=phases, weights=weights)
elif convert_tec:
    solsetout.makeSoltab('phase', axesNames=['pol', 'dir', 'ant', 'freq', 'time'], axesVals=[['XX'], directions, antennas, ax_freq, ax_time], vals=phases, weights=weights)
h5out.close()
