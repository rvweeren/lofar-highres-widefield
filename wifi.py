#!/usr/bin/env python
import configparser
import glob
import logging
import os
import subprocess
import sys
import traceback

import bdsf


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
try:
    # See if there is a running directory defined.
    os.chdir(os.path.expandvars("$RUNDIR"))
except FileNotFoundError:
    # If not, that's fine, run in the current directory.
    pass
CWD = os.getcwd()

# Set up logging stuff.
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()
logger.name = sys.argv[0]

logging.addLevelName(logging.INFO, "\033[1;32m%s\033[1;0m" % logging.getLevelName(logging.INFO))
logging.addLevelName(logging.WARNING, "\033[1;31m%s\033[1;0m" % logging.getLevelName(logging.WARNING))
logging.addLevelName(logging.ERROR, "\033[1;41m%s\033[1;0m" % logging.getLevelName(logging.ERROR))

# Read in the configuration file.
config = configparser.ConfigParser()
config.read(sys.argv[1])

logger.info('Checking inputs.')
if not os.path.isdir(config['data']['highres_data']):
    # Used for msoaicing at 0.2-0.3''.
    # die('High resolution data not found!')
    pass
if not os.path.isfile(config['solutions']['infield_sols_p']):
    die('Infield calibrator phase solutions not found!')
if not os.path.isfile(config['solutions']['infield_sols_ap']):
    die('Infield calibrator amp+phase solutions not found!')

mses = []
with open(config['data']['mslist']) as f:
    for l in f.readlines():
        mses.append(l.strip())

if config['data'].getboolean('do_apply_kms'):
    # This will apply the DIS2 solutions from the ddf-pipeline to arrive at DATA_DI_CORRECTED.
    if not os.path.exists(config['solutions']['kms_solsdir']):
        die('killMS solution directory not found!')
    logger.info('Converting kMS solutions to H5Parm.')
    for ms in mses:
        sols_npz = config['solutions']['kms_solsdir'] + '/' + ms + '/killMS.DIS2_full.sols.npz'
        sols_h5 = os.getcwd() + '/' + ms + '_DIS2_full.sols.h5'
        try:
            cmd = 'killMS2H5parm.py {h5:s} {npz:s} --nofulljones'.format(h5=sols_h5, npz=sols_npz)
            logger.info(cmd)
            subprocess.call(cmd, shell=True)
            cmd2 = 'addIS_to_h5.py {h5:s} {ms:s} --solset_in sol000 --solset_out sol001 --do_int_stations'.format(h5=sols_h5, ms=ms)
            logger.info(cmd2)
            subprocess.call(cmd2, shell=True)
        except Exception as e:
            traceback.print_exc()
            die()

    logger.info('Applying kMS solutions to MS.')
    for ms in mses:
        with open('apply_kms.parset', 'w') as f:
            sols = os.getcwd() + '/' + ms + '_DIS2_full.sols.h5'
            f.write('msin={ms:s} msin.datacolumn=DATA msout=. msout.datacolumn=DATA_DI_CORRECTED msout.storagemanager=dysco steps=[applykms] applykms.type=applycal applykms.steps=[p,a] applykms.parmdb={h5:s} applykms.solset={ss:s} applykms.p.correction=phase000 applykms.a.correction=amplitude000'.format(dc=config['data']['data_column'], ms=ms, h5=sols, ss='sol001').replace(' ', '\n'))
        cmd = 'DPPP apply_kms.parset'
        logger.info(cmd)
        subprocess.call(cmd, shell=True)

if config['data'].getboolean('do_subtract'):
    # Subtract sources outside a given region using the DDS3 solutions from the ddf-pipeline.
    # This is especially important with bright sources outside the FoV of the international stations,
    # but inside that of the Dutch stations.
    # Load the required settings.
    dc = config['subtract']['subtract_from']
    box = config['subtract']['boxfile']
    # Copy over the required files.
    path = config['subtract']['lotss_directory']
    logger.info('Copying over necessary LoTSS products.')
    import shutil
    from distutils.dir_util import copy_tree
    reqs = ['image_dirin_SSD_m.npy.ClusterCat.npy', 'DDS3_full_5038110493.005561_smoothed.npz', 'DDS3_full_slow_5038110493.005561_merged.npz', 'image_full_ampphase_di_m.NS.DicoModel', 'image_full_ampphase_di_m.NS.mask01.fits', 'SOLSDIR']
    for r in reqs:
        if os.path.isfile(path + '/' + r):
            shutil.copy2(path + '/' + r, os.getcwd() + '/')
        elif os.path.isdir(path + '/' + r):
            shutil.copytree(path + '/' + r, os.getcwd() + '/' + r)
    logger.info('Flagging international stations in all MS.')
    for ms in mses:
        cmd1 = 'backup_flagtable.py {:s}'.format(ms)
        cmd2 = 'DPPP flag_IS.parset msin={:s}'.format(ms)
        logger.info(cmd1)
        subprocess.call(cmd1, shell=True)
        logger.info(cmd2)
        subprocess.call(cmd2, shell=True)

    cmd = 'sub-sources-outside-region.py -b {:s} -m {:s} -c {:s} -f 1 -t 1 -p keepcenter'.format(box, config['data']['mslist'], dc)
    logger.info(cmd)
    subprocess.call(cmd, shell=True)

    logger.info('Restoring flags.')
    for ms in mses:
        cmd3 = 'restore_flagtable.py {:s}'.format(ms)
        logger.info(cmd3)
        subprocess.call(cmd3, shell=True)

if config['data'].getboolean('do_apply_infield'):
    if os.path.isfile(os.getcwd() + '/' + 'image_full_ampphase_di_m.NS_SUB.log'):
        logger.info('Applying infield calibrator solutions: DATA_SUB -> CORRECTED_DATA')
        dc = 'DATA_SUB'
    else:
        logger.info('Applying infield calibrator solutions: {:s} -> CORRECTED_DATA'.format(config['data']['data_column']))
        dc = config['data']['data_column']
    for ms in mses:
        sols_p = config['solutions']['infield_sols_p']
        if not config['solutions']['infield_phase_only']:
            sols_ap = config['solutions']['infield_sols_ap']
            try:
                with open('apply_infield_solutions.parset', 'w') as f:
                    f.write('msin={ms:s} msin.datacolumn={dc:s} msout=. msout.datacolumn=CORRECTED_DATA msout.storagemanager=dysco steps=[applyif1,applyif2] applyif1.type=applycal applyif1.parmdb={h51:s} applyif1.solset={ss1:s} applyif1.correction=phase000 applyif2.type=applycal applyif2.steps=[p,a] applyif2.parmdb={h52:s} applyif2.solset={ss2:s} applyif2.p.correction=phase000 applyif2.a.correction=amplitude000'.format(dc=dc, ms=ms, h51=sols_p, ss1=config['solutions']['infield_solset'], h52=sols_ap, ss2=config['solutions']['infield_solset']).replace(' ', '\n'))
                cmd = 'DPPP apply_infield_solutions.parset'
                logger.info(cmd)
                subprocess.call(cmd, shell=True)
            except Exception as e:
                traceback.print_exc()
                die()
        else:
            try:
                with open('apply_infield_solutions.parset', 'w') as f:
                    f.write('msin={ms:s} msin.datacolumn={dc:s} msout=. msout.datacolumn=CORRECTED_DATA msout.storagemanager=dysco steps=[applyif1] applyif1.type=applycal applyif1.parmdb={h51:s} applyif1.solset={ss1:s} applyif1.correction=phase000'.format(dc=dc, ms=ms, h51=sols_p, ss1=config['solutions']['infield_solset']).replace(' ', '\n'))
                os.system('DPPP apply_infield_solutions.parset')
            except Exception as e:
                traceback.print_exc()
                die()
else:
    logger.info('Infield solutions have been applied, skipping applycal step.')

tapered_images = glob.glob('wsclean_block*_taper*.fits')
if len(tapered_images):
    logger.info('Taper has already been created, skipping WSClean step.')
else:
    logger.info('Tapering data to target resolution of {:s}.'.format(config['image']['taper_full']))
    chan_out = (len(mses) // 4) + 1
    cmd = 'wsclean -j {:d} -mem {:d} -data-column {:s} -niter 0 -channels-out {:d} -weight briggs {:s} -size 1024 1024 -scale 0.35asec -make-psf -fit-psf -no-reorder -no-update-model-required -store-imaging-weights -taper-gaussian {:s}asec -name wsclean_taper *.{:s}'.format(int(config['image']['wsclean_ncpu']), int(config['image']['wsclean_mem']), config['image']['data_column'], chan_out, config['image']['robust_full'], config['image']['taper_full'], mses[0].split('.')[-1])
    logger.info(cmd)
    subprocess.call(cmd, shell=True)
    for i, ms in enumerate(mses):
        cmd = 'transfer_imaging_weight.py {:s}'.format(ms)
        logger.info(cmd)
        subprocess.call(cmd, shell=True)

if os.path.exists(os.getcwd() + '/image_dirin_SSD_init_natural.int.restored.fits'):
    logger.info('Initial widefield image already exists, not recreating.')
else:
    logger.info('Creating {:s}" widefield image.'.format(config['image']['taper_full']))
    cmd = 'DDF.py --Output-Name=image_dirin_SSD_init_natural --Data-MS={:s} --Deconv-PeakFactor 0.050000 --Data-ColName {ic:s} --Data-ChunkHours 4 --Parallel-NCPU=32 --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=10000 --Deconv-MaxMajorIter=3 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Mode Natural  --Image-NPix=25000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell {cell:f} --Facets-NFacets=7 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=. --Log-Memory 1 --GAClean-RMSFactorInitHMP 1.000000 --GAClean-MaxMinorIterInitHMP 10000.000000 --DDESolutions-SolsDir=SOLSDIR --Cache-Weight=reset --Output-Mode=Clean --Output-RestoringBeam 1.000000 --Weight-ColName="IMAGING_WEIGHT" --Freq-NBand=2 --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=10.00 --Selection-UVRangeKm=[5.0,2000.000000] --GAClean-MinSizeInit=10 --Output-RestoringBeam 1.000'.format(config['data']['mslist'], ic=config['image']['data_column'], cell=float(config['image']['cellsize_full']))
    logger.info(cmd)
    subprocess.call(cmd, shell=True)

if os.path.exists(os.getcwd() + '/image_dirin_SSD_init_natural.app.restored.fits.mask.fits'):
    logger.info('First mask already exists, not recreating.')
else:
    logger.info('Creating mask from initial image.')
    cmd = 'MakeMask.py --RestoredIm=image_dirin_SSD_init_natural.app.restored.fits --Th=7.5 --Box=50,2'
    logger.info(cmd)
    subprocess.call('MakeMask.py --RestoredIm=image_dirin_SSD_init_natural.app.restored.fits --Th=7.5 --Box=50,2', shell=True)

if os.path.exists(os.getcwd() + '/image_dirin_SSD_init_natural_m.int.restored.fits'):
    logger.info('Mask-cleaned image already exists, not recreating.')
else:
    logger.info('Cleaning deeper with mask.')
    cmd = 'DDF.py --Output-Name=image_dirin_SSD_init_natural_m --Data-MS={:s} --Deconv-PeakFactor 0.050000 --Data-ColName {ic:s} --Data-ChunkHours 4 --Parallel-NCPU=32 --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=10000 --Deconv-MaxMajorIter=5 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Mode Natural  --Image-NPix=25000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell {cell:f} --Facets-NFacets=7 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=. --Log-Memory 1 --GAClean-RMSFactorInitHMP 1.000000 --GAClean-MaxMinorIterInitHMP 10000.000000 --DDESolutions-SolsDir=SOLSDIR --Cache-Weight=reset --Output-Mode=Clean --Output-RestoringBeam 1.000000 --Weight-ColName="IMAGING_WEIGHT" --Freq-NBand=2 --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=10.00 --Selection-UVRangeKm=[5.0,2000.000000] --GAClean-MinSizeInit=10 --Output-RestoringBeam 1.000 --Mask-External=image_dirin_SSD.app.restored.fits.mask.fits --Predict-InitDicoModel=image_dirin_SSD_init_natural.DicoModel --Cache-Dirty=forceresidual'.format(config['data']['mslist'], ic=config['image']['data_column'], cell=float(config['image']['cellsize_full']))
    logger.info(cmd)
    subprocess.call(cmd, shell=True)

if os.path.exists(os.getcwd() + '/image_dirin_SSD_init_natural_m.app.restored.fits.mask.fits'):
    logger.info('Second mask already exists, not recreating.')
else:
    logger.info('Creating second mask.')
    cmd = 'MakeMask.py --RestoredIm=image_dirin_SSD_init_natural_m.app.restored.fits --Th=5 --Box=50,2'
    logger.info(cmd)
    subprocess.call(cmd, shell=True)

if os.path.exists(os.getcwd() + '/image_dirin_SSD_init_natural_m2.int.restored.fits'):
    logger.info('Final image already exists, not recreating.')
else:
    logger.info('Cleaning again with second mask.')
    cmd = 'DDF.py --Output-Name=image_dirin_SSD_init_natural_m2 --Data-MS={:s} --Deconv-PeakFactor 0.050000 --Data-ColName {ic:s} --Data-ChunkHours 4 --Parallel-NCPU=32 --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=10000 --Deconv-MaxMajorIter=3 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Mode Natural  --Image-NPix=25000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell {cell:f} --Facets-NFacets=7 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=. --Log-Memory 1 --GAClean-RMSFactorInitHMP 1.000000 --GAClean-MaxMinorIterInitHMP 10000.000000 --DDESolutions-SolsDir=SOLSDIR --Cache-Weight=reset --Output-Mode=Clean --Output-RestoringBeam 1.000000 --Weight-ColName="IMAGING_WEIGHT" --Freq-NBand=2 --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=10.00 --Selection-UVRangeKm=[5.0,2000.000000] --GAClean-MinSizeInit=10 --Output-RestoringBeam 1.000 --Mask-External=image_dirin_SSD_init_natural_m.app.restored.fits.mask.fits --Predict-InitDicoModel=image_dirin_SSD_init_natural_m.DicoModel --Cache-Dirty=forceresidual'.format(config['data']['mslist'], ic=config['image']['data_column'], cell=float(config['image']['cellsize_full']))
    logger.info(cmd)
    subprocess.call(cmd, shell=True)

logger.info('Making PyBDSF catalogue to select potential DDE calibrators.')
fitsname = 'image_dirin_SSD_init_natural_m2.int.restored.fits'
detectimage = 'image_dirin_SSD_init_natural_m2.app.restored.fits'
# Pull the reference frequency from the header.
fhead = fits.open('image_dirin_SSD_init_natural_m2.int.restored.fits')
restfrq = fhead[0].header['CRVAL3']
# Run PyBDSF with standard SKSP settings.
res = bdsf.process_image(fitsname, detection_image=detectimage, thresh_isl=4.0, thresh_pix=5.0, rms_box=(150, 15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60, 15), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=True, atrous_do=True, atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5, advanced_opts=True, blank_limit=None, frequency=restfrq)
# Write out a catalog.
res.write_catalog(outfile='skymodel_1asec_lbregion_pybdsf.csv', bbs_patches='source', catalog_type='gaul', format='csv')
res.write_catalog(outfile='skymodel_1asec_lbregion_pybdsf.bbs', bbs_patches='source', catalog_type='gaul', format='bbs')

logger.info('Pipeline finished successfully.')
sys.exit(0)

if config['mosaic'].getboolean('do_mosaic'):
    logger.info('0.2" mosaic requested. [NOT YET IMPLEMENTED]')
