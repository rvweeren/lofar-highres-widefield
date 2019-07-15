#!/usr/bin/env python
""" A widefield imaging pipeline for LOFAR HBA.
"""
import configparser
import glob
import logging
import os
import subprocess
import sys
import traceback

from astropy.io import fits

import bdsf
import numpy as np


def die(reason=''):
    ''' Stop the pipeline, reporting an error message if given.

    Args:
        reason (str): the reason for stopping the pipeline.
    Returns:
        None
    '''
    if reason:
        LOGGER.error(reason)
    else:
        LOGGER.error('Something went wrong for an unknown reason!')
    sys.exit(-1)


def get_beam(image):
    ''' Extracts the restoring beam from a FITS image and returns a list containing
    the major and minor axis, and the position angle.

    Args:
        image (str): image from which to extract the restoring beam.
    Returns:
        beam (list): a list with major axis, minor axis and position angle, all in degrees.
    '''
    f = fits.open(image)
    head = f[0].header
    bmaj = head['BMAJ']
    bmin = head['BMIN']
    bpa = head['BPA']
    temp = [bmaj, bmin, bpa]
    beam = list(np.round(temp, 10))

    return beam


def get_mslist():
    ''' Put all measurement sets from the mslist text file into a list.
    '''
    mses = []
    with open(CONFIG['data']['mslist']) as f:
        for l in f.readlines():
            mses.append(l.strip())
    return mses


def is_tapered():
    ''' Checks if the data has already been tapered.

    Returns:
        tapered (bool): `True` if tapering has been done, `False` otherwise.
    '''
    tapered_images = glob.glob('wsclean_taper*.fits')
    tapered = len(tapered_images) > 0
    return tapered


def run_pybdsf(fitsname, detectimage):
    ''' Run PyBDSF on an image, using standard SKSP settings.

    Two catalogues are written when finished:
        - A CSV catalogue with all columns present.
        - A BBS formatted catalogue, suitable for e.g. DPPP.

    Args:
        fitsname (str): path to the image from which to extract the fluxes.
        detectimage (str): path to the image on which to run source detection.
    Returns:
        None
    '''
    # Pull the reference frequency from the header.
    f = fits.open(fitsname)
    restfrq = f[0].header['CRVAL3']
    # Run PyBDSF with standard SKSP settings.
    res = bdsf.process_image(fitsname, detection_image=detectimage, thresh_isl=4.0, thresh_pix=5.0, rms_box=(150, 15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60, 15), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=True, atrous_do=True, atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5, advanced_opts=True, blank_limit=None, frequency=restfrq)
    # Write out a catalog.
    res.write_catalog(outfile='skymodel_1asec_lbregion_pybdsf.csv', bbs_patches='source', catalog_type='gaul', format='csv')
    res.write_catalog(outfile='skymodel_1asec_lbregion_pybdsf.bbs', bbs_patches='source', catalog_type='gaul', format='bbs')


# One should have run `genericpipeline.py -d -c pipeline.cfg LB-Split-Calibrators.parset` before running this.
# Two datasets must be present:
# - blocks of 10SB for the target field
# - a full bandwidth dataset for the infield calibrator
try:
    # See if there is a running directory defined.
    os.chdir(os.path.expandvars("$RUNDIR"))
except OSError:
    # If not, that's fine, run in the current directory.
    pass
CWD = os.getcwd()

# Set up logging stuff.
logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger()
LOGGER.name = sys.argv[0]

logging.addLevelName(logging.INFO, "\033[1;32m%s\033[1;0m" % logging.getLevelName(logging.INFO))
logging.addLevelName(logging.WARNING, "\033[1;31m%s\033[1;0m" % logging.getLevelName(logging.WARNING))
logging.addLevelName(logging.ERROR, "\033[1;41m%s\033[1;0m" % logging.getLevelName(logging.ERROR))

# Read in the configuration file.
CONFIG = configparser.ConfigParser()
CONFIG.read(sys.argv[1])

LOGGER.info('Checking inputs.')
if not os.path.isdir(CONFIG['data']['highres_data']):
    # Used for msoaicing at 0.2-0.3''.
    # die('High resolution data not found!')
    pass
if not os.path.isfile(CONFIG['solutions']['infield_sols_p']):
    die('Infield calibrator phase solutions not found!')
if not os.path.isfile(CONFIG['solutions']['infield_sols_ap']):
    die('Infield calibrator amp+phase solutions not found!')

MSES = get_mslist()

if CONFIG['data'].getboolean('do_apply_kms'):
    # This will apply the DIS2 solutions from the ddf-pipeline to arrive at DATA_DI_CORRECTED.
    if not os.path.exists(CONFIG['solutions']['kms_solsdir']):
        die('killMS solution directory not found!')
    LOGGER.info('Converting kMS solutions to H5Parm.')
    for ms in MSES:
        sols_npz = CONFIG['solutions']['kms_solsdir'] + '/' + ms + '/killMS.DIS2_full.sols.npz'
        sols_h5 = os.getcwd() + '/' + ms + '_DIS2_full.sols.h5'
        try:
            CMD = 'killMS2H5parm.py {h5:s} {npz:s} --nofulljones'.format(h5=sols_h5, npz=sols_npz)
            LOGGER.info(CMD)
            subprocess.call(CMD, shell=True)
            CMD2 = 'addIS_to_h5.py {h5:s} {ms:s} --solset_in sol000 --solset_out sol001 --do_int_stations'.format(h5=sols_h5, ms=ms)
            LOGGER.info(CMD2)
            subprocess.call(CMD2, shell=True)
        except Exception:
            traceback.print_exc()
            die()

    LOGGER.info('Applying kMS solutions to MS.')
    for ms in MSES:
        with open('apply_kms.parset', 'w') as f:
            sols = os.getcwd() + '/' + ms + '_DIS2_full.sols.h5'
            f.write('msin={ms:s} msin.datacolumn=DATA msout=. msout.datacolumn=DATA_DI_CORRECTED msout.storagemanager=dysco steps=[applykms] applykms.type=applycal applykms.steps=[p,a] applykms.parmdb={h5:s} applykms.solset={ss:s} applykms.p.correction=phase000 applykms.a.correction=amplitude000'.format(dc=CONFIG['data']['data_column'], ms=ms, h5=sols, ss='sol001').replace(' ', '\n'))
        CMD = 'DPPP apply_kms.parset'
        LOGGER.info(CMD)
        subprocess.call(CMD, shell=True)

if CONFIG['data'].getboolean('do_subtract'):
    # Subtract sources outside a given region using the DDS3 solutions from the ddf-pipeline.
    # This is especially important with bright sources outside the FoV of the international stations,
    # but inside that of the Dutch stations.
    # Load the required settings.
    dc = CONFIG['subtract']['subtract_from']
    box = CONFIG['subtract']['boxfile']
    # Copy over the required files.
    path = CONFIG['subtract']['lotss_directory']
    LOGGER.info('Copying over necessary LoTSS products.')
    import shutil
    reqs = ['image_dirin_SSD_m.npy.ClusterCat.npy', 'DDS3_full_5038110493.005561_smoothed.npz', 'DDS3_full_slow_5038110493.005561_merged.npz', 'image_full_ampphase_di_m.NS.DicoModel', 'image_full_ampphase_di_m.NS.mask01.fits', 'SOLSDIR']
    for r in reqs:
        if os.path.isfile(path + '/' + r):
            shutil.copy2(path + '/' + r, os.getcwd() + '/')
        elif os.path.isdir(path + '/' + r):
            shutil.copytree(path + '/' + r, os.getcwd() + '/' + r)
    LOGGER.info('Flagging international stations in all MS.')
    for ms in MSES:
        CMD1 = 'backup_flagtable.py {:s}'.format(ms)
        CMD2 = 'DPPP flag_IS.parset msin={:s}'.format(ms)
        LOGGER.info(CMD1)
        subprocess.call(CMD1, shell=True)
        LOGGER.info(CMD2)
        subprocess.call(CMD2, shell=True)

    CMD = 'sub-sources-outside-region.py -b {:s} -m {:s} -c {:s} -f 1 -t 1 -p keepcenter'.format(box, CONFIG['data']['mslist'], dc)
    LOGGER.info(CMD)
    subprocess.call(CMD, shell=True)

    LOGGER.info('Restoring flags.')
    for ms in MSES:
        CMD3 = 'restore_flagtable.py {:s}'.format(ms)
        LOGGER.info(CMD3)
        subprocess.call(CMD3, shell=True)

if CONFIG['control']['exitafter'] == 'subtract':
    LOGGER.info('Pipeline finished successfully.')
    sys.exit(0)

if CONFIG['data'].getboolean('do_apply_infield'):
    if os.path.isfile(os.getcwd() + '/' + 'image_full_ampphase_di_m.NS_SUB.log'):
        LOGGER.info('Applying infield calibrator solutions: DATA_SUB -> CORRECTED_DATA')
        dc = 'DATA_SUB'
    else:
        LOGGER.info('Applying infield calibrator solutions: {:s} -> CORRECTED_DATA'.format(CONFIG['data']['data_column']))
        dc = CONFIG['data']['data_column']
    for ms in MSES:
        sols_p = CONFIG['solutions']['infield_sols_p']
        if not CONFIG['solutions']['infield_phase_only']:
            sols_ap = CONFIG['solutions']['infield_sols_ap']
            try:
                with open('apply_infield_solutions.parset', 'w') as f:
                    f.write('msin={ms:s} msin.datacolumn={dc:s} msout=. msout.datacolumn=CORRECTED_DATA msout.storagemanager=dysco steps=[applyif1,applyif2] applyif1.type=applycal applyif1.parmdb={h51:s} applyif1.solset={ss1:s} applyif1.correction=phase000 applyif2.type=applycal applyif2.steps=[p,a] applyif2.parmdb={h52:s} applyif2.solset={ss2:s} applyif2.p.correction=phase000 applyif2.a.correction=amplitude000'.format(dc=dc, ms=ms, h51=sols_p, ss1=CONFIG['solutions']['infield_solset'], h52=sols_ap, ss2=CONFIG['solutions']['infield_solset']).replace(' ', '\n'))
                CMD = 'DPPP apply_infield_solutions.parset'
                LOGGER.info(CMD)
                subprocess.call(CMD, shell=True)
            except Exception:
                traceback.print_exc()
                die()
        else:
            try:
                with open('apply_infield_solutions.parset', 'w') as f:
                    f.write('msin={ms:s} msin.datacolumn={dc:s} msout=. msout.datacolumn=CORRECTED_DATA msout.storagemanager=dysco steps=[applyif1] applyif1.type=applycal applyif1.parmdb={h51:s} applyif1.solset={ss1:s} applyif1.correction=phase000'.format(dc=dc, ms=ms, h51=sols_p, ss1=CONFIG['solutions']['infield_solset']).replace(' ', '\n'))
                os.system('DPPP apply_infield_solutions.parset')
            except Exception:
                traceback.print_exc()
                die()
else:
    LOGGER.info('Infield solutions have been applied, skipping applycal step.')

if is_tapered():
    LOGGER.info('Taper has already been created, skipping WSClean step.')
else:
    LOGGER.info('Tapering data to target resolution of {:s}.'.format(CONFIG['image']['taper_full']))
    chan_out = (len(MSES) // 4) + 1
    # Having a small pixel scale is important here, such that all baselines are considered during
    # weighting. WSClean automatically ignores all baselines that provide resolutions higher than
    # the given pixel size.
    CMD = 'wsclean -j {:d} -mem {:d} -data-column {:s} -niter 0 -channels-out {:d} -weight briggs {:s} -size 1024 1024 -scale 0.05asec -minuvw-m 5000 -make-psf -fit-beam -no-reorder -no-update-model-required -store-imaging-weights -taper-gaussian {:s}asec -name wsclean_taper *.{:s}'.format(int(CONFIG['image']['wsclean_ncpu']), int(CONFIG['image']['wsclean_mem']), CONFIG['image']['data_column'], chan_out, CONFIG['image']['robust_full'], CONFIG['image']['taper_full'], MSES[0].split('.')[-1])
    LOGGER.info(CMD)
    subprocess.call(CMD, shell=True)
    for i, ms in enumerate(MSES):
        CMD = 'transfer_imaging_weight.py {:s}'.format(ms)
        LOGGER.info(CMD)
        subprocess.call(CMD, shell=True)
    #BEAM = get_beam('wsclean_taper-MFS-psf.fits')
    #DDF_RESTORING_BEAM = '[{maj:f},{min:f},{pa:f}'.format(BEAM[0], BEAM[1], BEAM[2])
    DDF_RESTORING_BEAM = '1.0'

if os.path.exists(os.getcwd() + '/image_dirin_SSD_init_natural.int.restored.fits'):
    LOGGER.info('Initial widefield image already exists, not recreating.')
else:
    LOGGER.info('Creating {:s}" widefield image.'.format(CONFIG['image']['taper_full']))
    CMD = 'DDF.py --Output-Name=image_dirin_SSD_init_natural --Data-MS={:s} --Deconv-PeakFactor 0.050000 --Data-ColName {ic:s} --Data-ChunkHours 4 --Parallel-NCPU=32 --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=10000 --Deconv-MaxMajorIter=3 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Mode Natural  --Image-NPix=25000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell {cell:f} --Facets-NFacets=7 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=. --Log-Memory 1 --GAClean-RMSFactorInitHMP 1.000000 --GAClean-MaxMinorIterInitHMP 10000.000000 --DDESolutions-SolsDir=SOLSDIR --Cache-Weight=reset --Output-Mode=Clean --Output-RestoringBeam {beam:s} --Weight-ColName="IMAGING_WEIGHT" --Freq-NBand=2 --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=10.00 --Selection-UVRangeKm=[5.0,2000.000000] --GAClean-MinSizeInit=10'.format(CONFIG['data']['mslist'], ic=CONFIG['image']['data_column'], cell=float(CONFIG['image']['cellsize_full'], beam=DDF_RESTORING_BEAM))
    LOGGER.info(CMD)
    subprocess.call(CMD, shell=True)

if os.path.exists(os.getcwd() + '/image_dirin_SSD_init_natural.app.restored.fits.mask.fits'):
    LOGGER.info('First mask already exists, not recreating.')
else:
    LOGGER.info('Creating mask from initial image.')
    CMD = 'MakeMask.py --RestoredIm=image_dirin_SSD_init_natural.app.restored.fits --Th=7.5 --Box=50,2'
    LOGGER.info(CMD)
    subprocess.call('MakeMask.py --RestoredIm=image_dirin_SSD_init_natural.app.restored.fits --Th=7.5 --Box=50,2', shell=True)

if os.path.exists(os.getcwd() + '/image_dirin_SSD_init_natural_m.int.restored.fits'):
    LOGGER.info('Mask-cleaned image already exists, not recreating.')
else:
    LOGGER.info('Cleaning deeper with mask.')
    CMD = 'DDF.py --Output-Name=image_dirin_SSD_init_natural_m --Data-MS={:s} --Deconv-PeakFactor 0.050000 --Data-ColName {ic:s} --Data-ChunkHours 4 --Parallel-NCPU=32 --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=10000 --Deconv-MaxMajorIter=5 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Mode Natural  --Image-NPix=25000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell {cell:f} --Facets-NFacets=7 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=. --Log-Memory 1 --GAClean-RMSFactorInitHMP 1.000000 --GAClean-MaxMinorIterInitHMP 10000.000000 --DDESolutions-SolsDir=SOLSDIR --Cache-Weight=reset --Output-Mode=Clean --Output-RestoringBeam {beam:s} --Weight-ColName="IMAGING_WEIGHT" --Freq-NBand=2 --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=10.00 --Selection-UVRangeKm=[5.0,2000.000000] --GAClean-MinSizeInit=10 --Mask-External=image_dirin_SSD.app.restored.fits.mask.fits --Predict-InitDicoModel=image_dirin_SSD_init_natural.DicoModel --Cache-Dirty=forceresidual'.format(CONFIG['data']['mslist'], ic=CONFIG['image']['data_column'], cell=float(CONFIG['image']['cellsize_full'], beam=DDF_RESTORING_BEAM))
    LOGGER.info(CMD)
    subprocess.call(CMD, shell=True)

if os.path.exists(os.getcwd() + '/image_dirin_SSD_init_natural_m.app.restored.fits.mask.fits'):
    LOGGER.info('Second mask already exists, not recreating.')
else:
    LOGGER.info('Creating second mask.')
    CMD = 'MakeMask.py --RestoredIm=image_dirin_SSD_init_natural_m.app.restored.fits --Th=5 --Box=50,2'
    LOGGER.info(CMD)
    subprocess.call(CMD, shell=True)

if os.path.exists(os.getcwd() + '/image_dirin_SSD_init_natural_m2.int.restored.fits'):
    LOGGER.info('Final image already exists, not recreating.')
else:
    LOGGER.info('Cleaning again with second mask.')
    CMD = 'DDF.py --Output-Name=image_dirin_SSD_init_natural_m2 --Data-MS={:s} --Deconv-PeakFactor 0.050000 --Data-ColName {ic:s} --Data-ChunkHours 4 --Parallel-NCPU=32 --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=10000 --Deconv-MaxMajorIter=3 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Mode Natural  --Image-NPix=25000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell {cell:f} --Facets-NFacets=7 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=. --Log-Memory 1 --GAClean-RMSFactorInitHMP 1.000000 --GAClean-MaxMinorIterInitHMP 10000.000000 --DDESolutions-SolsDir=SOLSDIR --Cache-Weight=reset --Output-Mode=Clean --Output-RestoringBeam {beam:s} --Weight-ColName="IMAGING_WEIGHT" --Freq-NBand=2 --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=10.00 --Selection-UVRangeKm=[5.0,2000.000000] --GAClean-MinSizeInit=10 --Mask-External=image_dirin_SSD_init_natural_m.app.restored.fits.mask.fits --Predict-InitDicoModel=image_dirin_SSD_init_natural_m.DicoModel --Cache-Dirty=forceresidual'.format(CONFIG['data']['mslist'], ic=CONFIG['image']['data_column'], cell=float(CONFIG['image']['cellsize_full'], beam=DDF_RESTORING_BEAM))
    LOGGER.info(CMD)
    subprocess.call(CMD, shell=True)

LOGGER.info('Making PyBDSF catalogue to select potential DDE calibrators.')

run_pybdsf(fitsname='image_dirin_SSD_init_natural_m2.int.restored.fits', detectimage='image_dirin_SSD_init_natural_m2.app.restored.fits')

LOGGER.info('Pipeline finished successfully.')
sys.exit(0)

if CONFIG['mosaic'].getboolean('do_mosaic'):
    LOGGER.info('0.2" mosaic requested. [NOT YET IMPLEMENTED]')
