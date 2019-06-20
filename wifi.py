#!/usr/bin/env python
import configparser
import glob
import logging
import os
import sys
import traceback

def die(reason=''):
    if reason:
        logger.error(reason)
    else:
        logger.error('Something went wrong for an unknown reason!')
    sys.exit(-1)

# One should have run `genericpipeline.py -d -c pipeline.cfg LB-Split-Calibrators.parset` before running this.
# Two datasets must be present:
# - blocks of 10SB for the target field
# - a full bandwidth dataset for the infield calibrator

CWD = os.getcwd()

# Set up logging stuff.
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()
logger.name = sys.argv[0]

logging.addLevelName( logging.INFO, "\033[1;32m%s\033[1;0m" % logging.getLevelName(logging.INFO))
logging.addLevelName( logging.WARNING, "\033[1;31m%s\033[1;0m" % logging.getLevelName(logging.WARNING))
logging.addLevelName( logging.ERROR, "\033[1;41m%s\033[1;0m" % logging.getLevelName(logging.ERROR))

# Read in the configuration file.
config = configparser.ConfigParser()
config.read(sys.argv[1])

logger.info('Checking inputs.')
if not os.path.isdir(config['data']['highres_data']):
    pass#die('High resolution data not found!')
if not os.path.isdir(config['solutions']['kms_solsdir']):
    pass#die('killMS DIS2 solutions not found!')
if not os.path.isfile(config['solutions']['infield_sols_p']):
    die('Infield calibrator phase solutions not found!')
if not os.path.isfile(config['solutions']['infield_sols_ap']):
    die('Infield calibrator amp+phase solutions not found!')

mses = []
with open(config['data']['mslist']) as f:
    for l in f.readlines():
        mses.append(l.strip())
'''
logger.info('Converting kMS solutions to H5Parm.')
for ms in mses:
    sols_npz = config['solutions']['kms_solsdir'] + '/' + ms + '/killMS.DIS2_full.sols.npz'
    sols_h5 = config['solutions']['kms_solsdir'] + '/' + ms + '/killMS.DIS2_full.sols.h5'
    try:
        print('killMS2H5parm.py {h5:s} {npz:s}'.format(h5=sols_h5, npz=sols_npz))
    except Exception as e:
        traceback.print_exc()
        die()
'''

if config['data'].getboolean('do_apply_infield'):
    logger.info('Applying infield calibrator solutions: DATA -> CORRECTED_DATA')
    for ms in mses:
        sols_p = config['solutions']['infield_sols_p']
        sols_ap = config['solutions']['infield_sols_ap']
        try:
            with open('apply_infield_solutions.parset', 'w') as f:
                f.write('msin={ms:s} msin.datacolumn={dc:s} msout=. msout.datacolumn=CORRECTED_DATA msout.storagemanager=dysco steps=[applyif1,applyif2] applyif1.type=applycal applyif1.parmdb={h51:s} applyif1.solset={ss1:s} applyif1.correction=phase000 applyif2.type=applycal applyif2.steps=[p,a] applyif2.parmdb={h52:s} applyif2.solset={ss2:s} applyif2.p.correction=phase000 applyif2.a.correction=amplitude000'.format(dc=config['data']['data_column'], ms=ms, h51=sols_p, ss1=config['solutions']['infield_solset'], h52=sols_ap, ss2=config['solutions']['infield_solset']).replace(' ', '\n'))
            os.system('DPPP apply_infield_solutions.parset')
        except Exception as e:
            traceback.print_exc()
            die()
else:
    logger.info('Infield solutions have been applied, skipping applycal step.')

logger.info('Tapering data to target resolution of {:s}.'.format(config['image']['taper_full']))
chan_out = (len(mses) // 4) + 1
print('wsclean -j {:d} -mem {:d} -data-column {:s} -niter 0 -weight briggs {:s} -size 1024 1024 -scale 0.025asec -store-imaging-weights -no-reorder -no-update-model-required -taper-gaussian {:s}asec -channels-out {:d} -name wsclean_taper {:s}'.format(int(config['image']['wsclean_ncpu']), int(config['image']['wsclean_mem']),config['image']['data_column'], config['image']['robust_full'], config['image']['taper_full'], chan_out, ' '.join(mses)))
os.system('wsclean -j {:d} -mem {:d} -data-column {:s} -niter 0 -weight briggs {:s} -size 1024 1024 -scale 0.025asec -store-imaging-weights -no-reorder -no-update-model-required -taper-gaussian {:s}asec -channels-out {:d} -name wsclean_taper {:s}'.format(int(config['image']['wsclean_ncpu']), int(config['image']['wsclean_mem']),config['image']['data_column'], config['image']['robust_full'], config['image']['taper_full'], chan_out, ' '.join(mses)))
for ms in mses:
    print('transfer_imgaging_weight.py {:s}'.format(ms))

#logger.info('Determining image weighting for lowest noise.')
#for r in range(-1.0, 1.0):
#    print('wsclean -no-update-model-required -minuv-l 80.0 -size 1024 1024 -reorder -weight briggs 0.6 -weighting-rank-filter 3 -clean-border 1 -mgain 0.8 -fit-beam -data-column {:s} -join-channels -channels-out 6 -padding 1.4 -auto-mask 2.5 -auto-threshold 1.0 -multiscale -multiscale-scales 0,4,8,16,32,64 -fit-spectral-pol 3 -pol i -baseline-averaging 0.266316109008 -name imagep1_b0.6 -scale 0.025arcsec -niter 0 {:s}'.format(config['data']['data_column'], mses[0]))

logger.info('Creating {:s}" widefield image.'.format(config['image']['taper_full']))
#print('DDF.py --Output-Name=image_dirin_SSD_init_natural --Data-MS=mslist.txt --Deconv-PeakFactor 0.050000 --Data-ColName {ic:s} --Data-ChunkHours 4 --Parallel-NCPU=32 --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=10000 --Deconv-MaxMajorIter=5 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Mode Natural  --Image-NPix=25000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell {cell:f} --Facets-NFacets=7 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=. --Log-Memory 1 --GAClean-RMSFactorInitHMP 1.000000 --GAClean-MaxMinorIterInitHMP 10000.000000 --DDESolutions-SolsDir=SOLSDIR --Cache-Weight=reset --Output-Mode=Clean --Output-RestoringBeam 1.000000 --Weight-ColName="IMAGING_WEIGHT" --Freq-NBand=2 --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=10.00 --Selection-UVRangeKm=[5.0,2000.000000] --GAClean-MinSizeInit=10 --Output-RestoringBeam 1.000'.format(ic=config['image']['data_column'], cell=float(config['image']['cellsize_full'])))
os.system('DDF.py --Output-Name=image_dirin_SSD_init_natural --Data-MS=mslist.txt --Deconv-PeakFactor 0.050000 --Data-ColName {ic:s} --Data-ChunkHours 4 --Parallel-NCPU=32 --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=10000 --Deconv-MaxMajorIter=5 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Mode Natural  --Image-NPix=25000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell {cell:f} --Facets-NFacets=7 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=. --Log-Memory 1 --GAClean-RMSFactorInitHMP 1.000000 --GAClean-MaxMinorIterInitHMP 10000.000000 --DDESolutions-SolsDir=SOLSDIR --Cache-Weight=reset --Output-Mode=Clean --Output-RestoringBeam 1.000000 --Weight-ColName="IMAGING_WEIGHT" --Freq-NBand=2 --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=10.00 --Selection-UVRangeKm=[5.0,2000.000000] --GAClean-MinSizeInit=10 --Output-RestoringBeam 1.000'.format(ic=config['image']['data_column'], cell=float(config['image']['cellsize_full'])))

logger.info('Pipeline finished successfully.')
sys.exit(0)
logger.info('Selecting potential DD calibration sources.')
import bdsf
fitsname = 'image_dirin_SSD_init_natural.int.restored.fits'
detectimage = 'image_dirin_SSD_init_natural.app.restored.fits'
# Pull the reference frequency from the header.
fhead = fits.open('image_dirin_SSD_init_natural.int.restored.fits')
restfrq= fhead[0].header['CRVAL3']
# Run PyBDSF with standard SKSP settings.
res = bdsf.process_image(fitsname, detection_image=detectimage, thresh_isl=4.0, thresh_pix=5.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=True, atrous_do=True, atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5, advanced_opts=True, blank_limit=None, frequency=restfrq)
# Write out a catalog.
res.write_catalog(outfile='skymodel_1asec_lbregion_pybdsf.csv', bbs_patches='source', catalog_type='gaul', format='csv')

if config['mosaic'].getboolean('do_mosaic'):
    logger.info('0.2" mosaic requested. [NOT YET IMPLEMENTED]')
