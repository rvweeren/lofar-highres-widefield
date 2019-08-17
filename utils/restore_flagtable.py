#!/usr/bin/env python
import sys

import casacore.tables as ct
import numpy as np

MS = sys.argv[1]

try:
    flags = np.load(sys.argv[2])
except:
    flags = np.load(MS + '.flagtable.npy')
tab = ct.table(MS, readonly=False)
tab.putcol('FLAG', flags)
