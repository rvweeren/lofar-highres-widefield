#!/usr/bin/env python
import sys

import casacore.tables as ct

ms = sys.argv[1]
try:
    tab = ct.taql('SELECT WEIGHT_SPECTRUM,IMAGING_WEIGHT_SPECTRUM FROM $ms')
except:
    print('Run wsclean with -store-imaging-weights on {:s} first!'.format(ms))
    sys.exit(-1)
# Add imaging columns, e.g. IMAGING_WEIGHT if they don't exist.
print('Updating IMAGING_WEIGHT of {:s}'.format(ms))
ct.addImagingColumns(ms)
try:
    tab = ct.taql('SELECT WEIGHT_SPECTRUM,IMAGING_WEIGHT_SPECTRUM,IMAGING_WEIGHT FROM $ms')
except:
    print('Run wsclean with -store-imaging-weights on {:s} first!'.format(ms))
    sys.exit(-1)
ws = tab.getcol('WEIGHT_SPECTRUM')
iws = tab.getcol('IMAGING_WEIGHT_SPECTRUM')
iw = (ws * iws)[:, :, 0]
print('Shape WEIGHT_SPECTRUM:')
print(ws.shape)
print('Shape IMAGING_WEIGHT_SPECTRUM:')
print(iws.shape)
print('Shape IMAGING_WEIGHT:')
print(iw.shape)

tab.putcol('IMAGING_WEIGHT', iw)
