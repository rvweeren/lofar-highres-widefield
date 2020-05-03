import argparse
import glob
import sys
import os
import re

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

def natural_sort( l ):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


def get_new_axes( mysoltab ):
    myvals = mysoltab.getValues()[0]
    axes_new = ['ant', 'freq', 'time']
    axes_vals = [ mysoltab.ant, mysoltab.freq, mysoltab.time ]
    myaxnames = mysoltab.getAxesNames()
    if 'dir' in myaxnames:
        axes_new = ['dir'] + axes_new
        axes_vals = [mysoltab.dir] + axes_vals
    if 'pol' in myaxnames:
        axes_new = ['pol'] + axes_new
        axes_vals = [mysoltab.pol] + axes_vals
    return( axes_new, axes_vals )

def main( mspath='', mssuffix='', h5parms='', soltabs2merge='phase,tec,amplitude', solsetin='sol000', h5out_name='', append_to_solset=None ):

    ## solution tables to merge
    print( soltabs2merge )
    soltabs = soltabs2merge.split(',')
    ordered_idx = []
    if 'phase' in soltabs:
        phs_idx = [ xx for xx in np.arange(len(soltabs)) if soltabs[xx] == 'phase' ]
        ordered_idx += phs_idx
        if 'tec' in soltabs:
            tec_idx = [ xx for xx in np.arange(len(soltabs)) if soltabs[xx] == 'tec' ]
            ordered_idx += tec_idx
    index = [ xx for xx in np.arange(len(soltabs)) if xx not in ordered_idx ]
    ordered_idx += index
    soltabs = np.array(soltabs)[ordered_idx].tolist()

    ## if '000' not in the name, add it
    soltabs = [ ss+'000' for ss in soltabs ]

    ## get time and freq resolution from data
    mssearch = os.path.join( mspath, mssuffix )
    mslist = natural_sort( glob.glob( mssearch ) )
    ms_first = mslist[0]
    ms_last = mslist[-1]

    print( ms_first, ms_last )

    print('Determining frequency grid...')
    ff = ct.taql('SELECT CHAN_FREQ, CHAN_WIDTH FROM ' + ms_first + '::SPECTRAL_WINDOW')
    freq_first = ff.getcol('CHAN_FREQ')[0][0]
    freq_spacing = ff.getcol('CHAN_WIDTH')[0][0]
    ff.close()

    fl = ct.taql('SELECT CHAN_FREQ, CHAN_WIDTH FROM ' + ms_last + '::SPECTRAL_WINDOW')
    freq_last = fl.getcol('CHAN_FREQ')[0][0]
    print(freq_first, freq_last, freq_spacing)
    ax_freq = np.arange(freq_first, freq_last + freq_spacing, freq_spacing)
    print( len(ax_freq) )

    print('Determining time grid...')
    ## to keep size down, determine from first h5
    if os.path.isfile(h5out_name):
        h5 = h5parm.h5parm(h5out_name)
    else:
        h5 = h5parm.h5parm(h5parms[0])
    ss = h5.getSolset(solsetin)
    soltab_names = ss.getSoltabNames()
    mint = []
    maxt = []
    tint = []
    for myst in soltab_names:
        tmp = ss.getSoltab(myst)
        times = tmp.time
        mint.append(np.min(times))
        maxt.append(np.max(times))
        tint.append(times[1]-times[0])
    smallest_idx = np.where(tint == np.min(tint))[0][0]
    ax_time = np.arange( mint[smallest_idx], maxt[smallest_idx]+tint[smallest_idx], tint[smallest_idx] )
    h5.close()

#    tf = ct.taql('SELECT TIME FROM ' + ms_first )
#    times = tf.getcol('TIME')
#    tf.close()
#    ax_time = np.unique( times )
    print( ax_time[0], ax_time[-1], np.float(ax_time[1])-np.float(ax_time[0]) )
    print( len( ax_time ) )

    ## get some information from the first h5parm
    h5 = h5parm.h5parm(h5parms[0])
    ss = h5.getSolset(solsetin)
    st = ss.getSoltab(soltabs[0])
    vals = st.getValues()[0]
    polarizations = st.pol
    antennas = st.ant
    ## use measurement set freq, times
    phases = np.zeros((len(polarizations), 1, len(antennas), len(ax_freq), len(ax_time)))
    phs_axes, phs_axes_vals = get_new_axes( st )
    if 'amplitude000' in soltabs:
        amps = np.ones((len(polarizations), 1, len(antennas), len(ax_freq), len(ax_time)))
        st = ss.getSoltab('amplitude000')
        vals= st.getValues()[0]
        amp_axes, amp_axes_vals = get_new_axes( st )

    ## antenna information
    antennas = st.getAxisValues('ant')
    ss_antennas = ss.obj.antenna.read()

    ## initialize the new h5 if it's not already there
    h5out = h5parm.h5parm(h5out_name, readonly=False)
    if append_to_solset is not None:
        solsetout = h5out.getSolset(append_to_solset)
    else:
        # Try to make a new one.
        solsetout = h5out.makeSolset('sol000')
    ## antenna information
    antennasout = solsetout.getAnt()
    antennatable = solsetout.obj._f_get_child('antenna')
    antennatable.append(ss.obj.antenna.read())

    h5.close()

    ## initialize directions and source list
    directions = []
    sourcelist = []

    for i, h5p in enumerate(h5parms):
        print('Processing direction for ' + h5p)
        ## open the h5parm and get source information
        h5 = h5parm.h5parm(h5p)
        ss = h5.getSolset(args.solsetin)
        source_coords = ss.getSou()['POINTING']
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

        ## start looping through soltabs
        for soltab in soltabs:
            st = ss.getSoltab(soltab)
            tmp_vals = st.getValues()[0]
            axes_new, axes_vals = get_new_axes( st )
            vals = reorderAxes( tmp_vals, st.getAxesNames(), axes_new )
            ## add phases to the total phase correction for this direction:
            if soltab == 'phase000':
                ## interpolate in time and frequency
                tmp_vals = interp_along_axis( vals, st.getAxisValues('time'), ax_time, -1 )
                vals = interp_along_axis( tmp_vals, st.getAxisValues('freq'), ax_freq, -2 )
                ## note that phase should have polarisation information
                if 'dir' in axes_new:
                    ## append if idx > len(phases)
                    if idx == 0:
                        phases[ :, idx, :, :, :] += vals[0, ...]
                    else:
                        phases = np.append(phases, vals, axis=1)
                else:
                    if idx == 0:
                        phases[ :, idx, :, :, :] += vals
                    else:
                        phases = np.append(phases, vals, axis=1)
            elif soltab == 'tec000':
                freqs = ax_freq.reshape(1,1,1,-1,1)
                tecphase = (-8.4479745e9 * vals / freqs)
                tp = interp_along_axis( tecphase, st.getAxisValues('time'), ax_time, -1)
                if 'pol' in phs_axes:
                    polidx = np.arange(len(polarizations))
                else:
                    polidx = 0
                if 'dir' in axes_new:
                    phases[ polidx, idx, :, :, :] += tp[0, ...]
                else:
                    phases[ polidx, idx, :, :, :] += tp
            elif soltab == 'amplitude000':
                tmp_vals = interp_along_axis( vals, st.getAxisValues('time'), ax_time, -1 )
                vals = interp_along_axis( tmp_vals, st.getAxisValues('freq'), ax_freq, -2 )
                if 'dir' in axes_new:
                    if idx == 0:
                        amps[ :, idx, :, :, :] += vals[0, ...]
                    else:
                        amps = np.append( amps, vals, axis=1 )
                else:
                    if idx == 0:
                        amps[ :, idx, :, :, :] *= vals
                    else:
                        amps = np.append( amps, vals, axis=1)
        h5.close()

    # Create the output h5parm.
    ## phases first
    weights = np.ones(phases.shape)
    sources = np.array(sourcelist, dtype=[('name', 'S128'), ('dir', '<f4', (2,))])
    solsetout.obj.source.append(sources)
    axes_vals = [polarizations, directions, antennas, ax_freq, ax_time]
    solsetout.makeSoltab('phase', axesNames=phs_axes, axesVals=axes_vals, vals=phases, weights=weights)
    ## and now amps
    weights = np.ones(amps.shape)
    solsetout.makeSoltab('amplitude', axesNames=amp_axes, axesVals=axes_vals, vals=amps, weights=weights)
    h5out.close()
    print( 'done.' )
    return


if __name__ == "__main__":

    parser =argparse.ArgumentParser()
    parser.add_argument('--mspath', dest='mspath', help='Path to the directory with easurement sets to pull frequency axis from, when converting TEC to phase.')
    parser.add_argument('--mssuffix', dest='mssuffix', default='ms', help='Suffix of your measurement sets, e.g. MS or ms.')
    parser.add_argument('--h5parms', dest='h5parms', nargs='+', help='Input H5parms to merge as directions, where each h5parm is one direction.')
    parser.add_argument('--soltabs', dest='soltabs2merge', help='SolTab of the H5parms to merge.', default='phase,tec,amplitude')
    parser.add_argument('--solset-in', dest='solsetin', help='SolSet to take the soltab from.', default='sol000')
    parser.add_argument('--h5parm-out', dest='h5out', help='Output H5parm with all directions present.')
    parser.add_argument('--append-to-solset', dest='append_to_solset', default=None, help='Append the new soltab to the given solset instead of creating a new one.')
    args = parser.parse_args()

    main(  mspath=args.mspath, mssuffix=args.mssuffix, h5parms=args.h5parms, soltabs2merge=args.soltabs2merge, solsetin=args.solsetin, h5out_name=args.h5out, append_to_solset=args.append_to_solset )
