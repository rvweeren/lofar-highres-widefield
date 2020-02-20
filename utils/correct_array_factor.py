#!/usr/bin/env python
import os, sys
import numpy as np
import losoto
import losoto.lib_operations
import glob
from astropy.io import fits
import pyrap.tables as pt
import os.path
from losoto import h5parm
import logging
import fnmatch

from lofar.stationresponse import stationresponse

def beamcor(ms):
    """  
    correct a ms for the beam in the phase center (array_factor only)
    """
    losoto = 'losoto'
    taql = 'taql'
    H5name = create_beamcortemplate(ms)
    parset = create_losoto_beamcorparset(ms)

    losotolofarbeam(H5name, 'phase000', ms, useElementResponse=False, useArrayFactor=True, useChanFreq=True)
    losotolofarbeam(H5name, 'amplitude000', ms, useElementResponse=False, useArrayFactor=True, useChanFreq=True)   

    fixbeam_ST001(H5name)

    cmdlosoto = losoto + ' ' + H5name + ' ' + parset
    print cmdlosoto
    os.system(cmdlosoto)
    
    cmd = 'DPPP numthreads=32 msin=' + ms + ' msin.datacolumn=DATA msout=. '
    cmd += 'msin.weightcolumn=WEIGHT_SPECTRUM '
    cmd += 'msout.datacolumn=CORRECTED_DATA steps=[ac1,ac2] msout.storagemanager=dysco '
    cmd += 'ac1.parmdb='+H5name + ' ac2.parmdb='+H5name + ' '
    cmd += 'ac1.type=applycal ac2.type=applycal '
    cmd += 'ac1.correction=phase000 ac2.correction=amplitude000 ac2.updateweights=True '
    print 'DPPP applycal:', cmd
    os.system(cmd)
   
    os.system(taql + " 'update " + ms + " set DATA=CORRECTED_DATA'")
    # Add beam correction keyword here.
    # This script only applies the array factor and assumes the element beam was corrected already.
    # Valid values are Element, ArrayFactor or Full.
    t = pt.table(ms, readonly=False)
    t.putcolkeywords('DATA', {'LOFAR_APPLIED_BEAM_MODE':'Full'})

    t2 = pt.table(ms + '::FIELD')
    phasedir = t2.getcol('PHASE_DIR').squeeze()
    t2.close()
    beamdir = t.getcolkeyword('DATA', 'LOFAR_APPLIED_BEAM_DIR')
    # Right ascension in radians is set in m0
    # Declination in radians is set in m1
    beamdir['m0']['value'] = phasedir[0]
    beamdir['m1']['value'] = phasedir[1]
    t.putcolkeywords('DATA', {'LOFAR_APPLIED_BEAM_DIR': beamdir})
    
    return

def create_beamcortemplate(ms):
  """  
  create a DPPP gain H5 template solutution file that can be filled with losoto
  """
  H5name = ms + '_templatejones.h5'   

  cmd = 'DPPP numthreads=32 msin=' + ms + ' msin.datacolumn=DATA msout=. '
  cmd += 'msin.modelcolumn=DATA '
  cmd += 'steps=[ddecal] ddecal.type=ddecal '
  cmd += 'ddecal.maxiter=1 ddecal.usemodelcolumn=True ddecal.nchan=1 '
  cmd += 'ddecal.mode=complexgain ddecal.h5parm=' + H5name  + ' '
  cmd += 'ddecal.solint=10'

  print cmd
  os.system(cmd)

  return H5name

def fixbeam_ST001(H5name):

   H5 = h5parm.h5parm(H5name, readonly=False)

   ants = H5.getSolset('sol000').getAnt().keys()
   antsrs = fnmatch.filter(ants,'RS*')

   if 'ST001' in ants:

     amps    = H5.getSolset('sol000').getSoltab('amplitude000').getValues()
     ampvals = H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0]
     phasevals = H5.getSolset('sol000').getSoltab('phase000').getValues()[0]

     idx = np.where(amps[1]['ant'] == 'ST001')[0][0]
     idxrs = np.where(amps[1]['ant'] == antsrs[0])[0][0]
     #idx106 = np.where(amps[1]['ant'] == 'RS106HBA')[0][0]
     #idx305 = np.where(amps[1]['ant'] == 'RS305HBA')[0][0]
     #idx508 = np.where(amps[1]['ant'] == 'RS508HBA')[0][0]
     #idx406 = np.where(amps[1]['ant'] == 'RS406HBA')[0][0]
     #idxnonecheck = np.where(amps[1]['ant'] == 'blaHBA')[0][0]

     #print idx205
     #print '$$$$$$$$$$$$$$$$$$$$$$$$'
     #print '$$$$$$$$$$$$$$$$$$$$$$$$'
     #print '$$$$$$$$$$$$$$$$$$$$$$$$'

     #ampvals[:,:, idx, 0,:] = 1.0  # set amplitude to 1.
     #phasevals[:,:, idx, 0,:] = 0.0 # set phase to 0.

     ampvals[:,:, idx, 0,:] = ampvals[:,:, idxrs, 0,:]
     phasevals[:,:, idx, 0,:] = 0.0

     H5.getSolset('sol000').getSoltab('amplitude000').setValues(ampvals)
     H5.getSolset('sol000').getSoltab('phase000').setValues(phasevals)

   H5.close()

   return

def create_losoto_beamcorparset(ms):
    """  
    Create a losoto parset to fill the beam correction values'.
    """
    parset = 'losotobeam.parset'
    os.system('rm -f ' + parset)
    f=open(parset, 'w') 

    f.write('pol = [XX,YY]\n')
    f.write('soltab = [sol000/*]\n\n\n')

    #f.write('[beam]\n')
    #f.write('ms = %s\n' % ms)
    #f.write('operation = LOFARBEAM\n')
    #f.write('useElementResponse = False\n')
    #f.write('useArrayFactor = True\n')
    #f.write('useChanFreq = True\n\n\n')

    f.write('[plotphase]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/phase000]\n')
    f.write('axesInPlot = [time,freq]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [-0.5,0.5]\n')
    f.write('prefix = plotlosoto%s/phases_beam\n' % ms)
    f.write('refAnt = CS001HBA0\n\n\n')

    f.write('[plotamp]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/amplitude000]\n')
    f.write('axesInPlot = [time,freq]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [0.2,1]\n')
    f.write('prefix = plotlosoto%s/amplitudes_beam\n' %ms)

    f.close()
    return parset

def losotolofarbeam(parmdb, soltabname, ms, inverse=False, useElementResponse=True, useArrayFactor=True, useChanFreq=True):
    """
    Do the beam correction via this imported losoto operation
    """

    H5 = h5parm.h5parm(parmdb, readonly=False)
    soltab = H5.getSolset('sol000').getSoltab(soltabname)

    #t = pt.table(ms)
    sr = stationresponse(ms, inverse, useElementResponse, useArrayFactor, useChanFreq)

    numants = pt.taql('select gcount(*) as numants from '+ms+'::ANTENNA').getcol('numants')[0]
    times = soltab.getAxisValues('time')

    for vals, coord, selection in soltab.getValuesIter(returnAxes=['ant','time','pol','freq'], weight=False):
        vals = losoto.lib_operations.reorderAxes( vals, soltab.getAxesNames(), ['ant','time','freq','pol'] )

        for stationnum in range(numants):
            logging.debug('Working on station number %i' % stationnum)
            for itime, time in enumerate(times):
                beam = sr.evaluateStation(time=time, station=stationnum)
                # Reshape from [nfreq, 2, 2] to [nfreq, 4]
                beam = beam.reshape(beam.shape[0], 4)

                if soltab.getAxisLen('pol') == 2:
                    beam = beam[:,[0,3]] # get only XX and YY

                if soltab.getType() == 'amplitude':
                    vals[stationnum, itime, :, :] = np.abs(beam)
                elif soltab.getType() == 'phase':
                    vals[stationnum, itime, :, :] = np.angle(beam)
                else:
                    logging.error('Beam prediction works only for amplitude/phase solution tables.')
                    return 1

        vals = losoto.lib_operations.reorderAxes( vals, ['ant','time','freq','pol'], [ax for ax in soltab.getAxesNames() if ax in ['ant','time','freq','pol']] )
        soltab.setValues(vals, selection)

    H5.close()
    return

if __name__ == '__main__':
    msin = sys.argv[1]
    beamcor(msin)
