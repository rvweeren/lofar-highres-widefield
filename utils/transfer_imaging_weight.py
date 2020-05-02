#!/usr/bin/env python
import argparse
import sys

import casacore.tables as ct


def main(ms):
    ''' Transfers a taper made with WSClean to IMAGING_WEIGTH such that it can be used in DDFacet.

    Args:
        ms (str): Measurement Set to operate on.

    Returns:
        None
    '''
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
    print('Reading WEIGHT_SPECTRUM and IMAGING_WEIGHT_SPECTRUM')
    ws = tab.getcol('WEIGHT_SPECTRUM')
    iws = tab.getcol('IMAGING_WEIGHT_SPECTRUM')
    iw = (ws * iws)[:, :, 0]
    print('Shape WEIGHT_SPECTRUM:')
    print(ws.shape)
    print('Shape IMAGING_WEIGHT_SPECTRUM:')
    print(iws.shape)
    print('Shape IMAGING_WEIGHT:')
    print(iw.shape)

    print('Writing taper to IMAGING_WEIGHT')
    tab.putcol('IMAGING_WEIGHT', iw)
    print('Taper written to IMAGING_WEIGHT')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('ms', help='Measurement Set to transfer weights in.')
    args = parser.parse_arguments()
    main(args.ms)
