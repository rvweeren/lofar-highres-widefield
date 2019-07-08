import sys

import casacore.tables as ct
import numpy as np

MS = sys.argv[1]

tab = ct.table(MS)
flags = tab.getcol('FLAG')
np.save(MS + '.flagtable', flags)
