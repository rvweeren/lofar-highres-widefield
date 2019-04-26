#!/usr/bin/env python
from __future__ import division
from astropy.io import fits

import argparse
import logging
import os
import subprocess
import time

import casacore.tables as ct


def do_dppp(parset):
    logging.info('Starting DPPP.')
    time_start = time.time()
    os.system('DPPP {parset}'.format(parset=parset))
    time_end = time.time()
    logging.info('DPPP took {} seconds.'.format(time_end - time_start))

def make_parset_ampphase(ms, solint=1, nchan=1, iteration=1, PARSET='selfcal_ampphase.parset'):
    if os.path.isfile(PARSET):
        os.remove(PARSET)
    PHASE_PARSET = '''msin = {msfile}
msin.datacolumn = CORRECTED_DATA_PHASE_ONLY
msout = .
msout.datacolumn = CORRECTED_DATA_AMPPHASE
msout.storagemanager = dysco

steps = [ampcal]

ampcal.type           = gaincal
ampcal.caltype        = diagonal
ampcal.usemodelcolumn = true
ampcal.parmdb         = {msfile}_selfcal_ap{i}.h5
ampcal.solint         = {solint}
ampcal.nchan          = {nchan}
ampcal.applysolution  = true
ampcal.uvmmin         = 40000'''
    with open(PARSET, 'w') as f:
        f.write(PHASE_PARSET.format(msfile=ms, solint=solint, nchan=nchan, i=iteration))
    return PARSET


def make_parset_losoto():
    pass


def make_parset_phaseonly(ms, solint=1, nchan=1, iteration=1, PARSET='selfcal_phaseonly.parset'):
    if os.path.isfile(PARSET):
        os.remove(PARSET)
    PHASE_PARSET = '''msin = {msfile}
msout = .
msout.datacolumn = CORRECTED_DATA_PHASE_ONLY
msout.storagemanager = dysco

steps = [phasecal]

phasecal.type           = gaincal
phasecal.caltype        = phaseonly
phasecal.usemodelcolumn = true
phasecal.parmdb         = {msfile}_selfcal_p{i}.h5
phasecal.solint         = {solint}
phasecal.nchan          = {nchan}
phasecal.applysolution  = true
phasecal.uvmmin         = 40000'''
    with open(PARSET, 'w') as f:
        f.write(PHASE_PARSET.format(msfile=ms, solint=solint, nchan=nchan, i=iteration))
    return PARSET


def run_wsclean(ms, colname='DATA', name='image', j=32, mem=75, scale='0.025asec', size=1024, baseline_averaging=100, no_update_model_required=True, niter=5000, auto_mask=5,\
    auto_threshold=1.5, local_rms=True, weight='briggs', robust=-1, channels_out=12, join_channels=True, stop_negative=True, minuvw_m=30000, taper_gaussian=0.0,\
    predict=False):
    time_start = time.time()
    command = 'wsclean -data-column {colname:s} -name {imgname:s} -scale {scale:s} -size {size:d} {size:d} -niter {niter:d} '.format(colname=colname, imgname=name, size=size, niter=niter, scale=scale)
    if not predict:
        if no_update_model_required:
            command += '-no-update-model-required '
        if join_channels:
            command += '-join-channels '
        if baseline_averaging:
            command += '-baseline-averaging {blavg:d} '.format(blavg=baseline_averaging)
        if stop_negative:
            command += '-stop-negative '
    if weight == 'briggs':
        command += '-weight briggs {robust:f} '.format(robust=robust) 
    if auto_mask:
        command += '-auto-mask {auto_mask:f} '.format(auto_mask=auto_mask)
    if auto_threshold:
        command += '-auto-threshold {autothresh:f} '.format(autothresh=auto_threshold)
    if local_rms:
        command += '-local-rms '
    if channels_out:
        command += '-channels-out {chano:d} '.format(chano=channels_out)
    if minuvw_m:
        command += '-minuvw-m {minuvwm:f} '.format(minuvwm=minuvw_m)
    if taper_gaussian:
        command += '-taper-gaussian {taper:f} '.format(taper=taper_gaussian)
    command += ms
    logging.info('Running: ' + command)
    os.system(command)
    time_end = time.time()
    logging.info('WSClean took {} seconds.'.format(time_end - time_start))

def measure_statistic(fitsname):
    # Assume a single image (no cube).
    img = fits.open(fitsname)[0].data.squeeze()
    imin = img.min()
    imax = img.max()
    statistic = abs(imax / imin)
    return statistic

def main(ms, solint_phase, nchan_phase, solint_amp, nchan_amp, image_name='image', colname='DATA'):
    if not os.path.isdir(ms):
        logging.warning('Measurement set not found!')
        return -1

    if (len(solint_phase) != len(nchan_phase)) or (len(solint_amp) != len(nchan_amp)):
        logging.warning('Solint and nchan lists should have the same length!')
    logging.basicConfig(level=logging.INFO)
    statistics = []
    if not os.path.isfile(image_name + '_initial-MFS-image.fits'):
        logging.info('MAKING INITIAL IMAGE')
        run_wsclean(ms, j=32, mem=75, scale='0.025asec', baseline_averaging=290, no_update_model_required=True, niter=5000, auto_mask=5, local_rms=True, auto_threshold=1.5, weight='briggs', robust=-1, channels_out=12, join_channels=True, stop_negative=True, size=1024, minuvw_m=30000, colname=colname, name=image_name+'_initial')
    else:
        logging.info('Initial image exists, skipping WSClean image step.')

    statistic_current = measure_statistic(image_name + '_initial-MFS-image.fits')
    statistics.append(statistic_current)
    if not os.path.isfile(ms+'_has_predict_initial'):
        logging.info('PREDICTING INITIAL IMAGE')
        run_wsclean(ms, j=32, mem=75, scale='0.025asec', baseline_averaging=290, no_update_model_required=True, niter=5000, auto_mask=5, local_rms=True, auto_threshold=1.5, weight='briggs', robust=-1, channels_out=12, join_channels=True, stop_negative=True, size=1024, minuvw_m=30000, colname=colname, name=image_name+'_initial', predict=True)
        subprocess.call('touch {ms}_has_predict_initial'.format(ms=ms), shell=True)
    else:
        logging.info('Initial model has been predicted, skipping WSClean predict step.')
        logging.info('Initalizing statistic to 1.')
        statistic_current = 1.

    # Start iterating phase-only self-calibration cycles.
    logging.info('Starting phase-only self-calibration.')
    for j, (ss, cc) in enumerate(zip(solint_phase[0], nchan_phase[0])):
        i = 0
        s = ss
        c = cc
        if (c == -1) or (s == -1):
            statistic_new = 1.
            logging.info('Negative solint or nchan, skipping phase-only!')
            break
        logging.info('Starting phase only solve with solint {:d} and nchan {:d}'.format(s, c))
        done = False
        while not done:
            logging.info('Iteration: ' + str(i))
            statistic_current = statistic_new
            if os.path.isfile('selfcal_p{iteration}.h5'.format(iteration=i)):
                # This iteration has been run already; skip to the next.
                logging.info('DPPP iteration {iteration} has been run already, skipping DPPP.'.format(iteration=i))
            else:
                # Solution file does not exist, run calibration.
                # Create parset and run DPPP.
                parset = make_parset_phaseonly(ms, solint=s, nchan=c, iteration=i)
                do_dppp(parset)

            # Image self-calibrated data, check if the image improved and predict back.
            if not os.path.isfile(image_name + '_p{iteration}-MFS-image.fits'.format(iteration=i)):
                run_wsclean(ms, j=32, mem=75, scale='0.025asec', baseline_averaging=290, no_update_model_required=True, niter=5000, auto_mask=5, local_rms=True, auto_threshold=1.5, weight='briggs', robust=-1, channels_out=12, join_channels=True, stop_negative=True, size=1024, minuvw_m=30000, colname='CORRECTED_DATA_PHASE_ONLY', name=image_name + '_p{iteration}'.format(iteration=i), predict=False)
            else:
                logging.info('Iteration p{iteration} has been imaged already, skipping WSClean.'.format(iteration=i))
            statistic_new = measure_statistic(image_name + '_p{iteration}-MFS-image.fits'.format(iteration=i))
            statistics.append(statistic_new)
            logging.info('All statistics: ' + ', '.join([str(stats) for stats in statistics]))
            logging.info('Current statistic: ' + str(statistic_current))
            logging.info('New statistic: ' + str(statistic_new))
            if (statistic_new / statistic_current) > 1.02:
                # Image quality improved, predict back.
                logging.info('Image quality improved.')
                stall = 0
            elif statistic_new < statistic_current:
                logging.info('Image quality did not improve, stopping phase-only iterations at this solint.')
                logging.info('Best iteration: {:d}'.format(i-1))
                done = True
            elif (statistic_new / statistic_current) < 1.02:
                logging.info('Image improved less than 2%.')
                stall += 1
                if stall > 3:
                    logging.info('Image improved less than 2% three successive times, selfcal is likely stalling.')
                    done = True
            logging.info('Predicting back to MODEL_DATA.')
            run_wsclean(ms, j=32, mem=75, scale='0.025asec', baseline_averaging=290, no_update_model_required=True, niter=5000, auto_mask=5, local_rms=True, auto_threshold=1.5, weight='briggs', robust=-1, channels_out=12, join_channels=True, stop_negative=True, size=1024, minuvw_m=30000, colname='CORRECTED_DATA_PHASE_ONLY', name=image_name + '_p{iteration}'.format(iteration=i), predict=True)
            subprocess.call('touch {ms}_has_predict_phase_p{iteration}'.format(ms=ms, iteration=i), shell=True)

    # Start iterating amplitude+phase self-calibration cycles.
    logging.info('Starting amplitude+phase self-calibration.')
    for j, (ss, cc) in enumerate(zip(solint_amp[0], nchan_amp[0])):
        i = 0
        s = ss
        c = cc
        if (c == -1) or (s == -1):
            statistic_new = 1.
            logging.info('Negative solint or nchan, skipping amplitude+phase!')
            break
        logging.info('Starting amplitude+phase solve with solint {:d} and nchan {:d}'.format(s, c))
        done = False
        while not done:
            logging.info('Iteration: ' + str(i))
            statistic_current = statistic_new
            if os.path.isfile('selfcal_ap{iteration}.h5'.format(iteration=i)):
                # This iteration has been run already; skip to the next.
                logging.info('DPPP iteration {iteration} has been run already, skipping DPPP.'.format(iteration=i))
            else:
                # Solution file does not exist, run calibration.
                # Create parset and run DPPP.
                parset = make_parset_ampphase(ms, solint=s, nchan=c, iteration=i)
                do_dppp(parset)

            # Image self-calibrated data, check if the image improved and predict back.
            if not os.path.isfile(image_name + '_ap{iteration}-MFS-image.fits'.format(iteration=i)):
                run_wsclean(ms, j=32, mem=75, scale='0.025asec', baseline_averaging=290, no_update_model_required=True, niter=5000, auto_mask=5, local_rms=True, auto_threshold=1.5, weight='briggs', robust=-1, channels_out=12, join_channels=True, stop_negative=True, size=1024, minuvw_m=30000, colname='CORRECTED_DATA_PHASE_ONLY', name=image_name + '_ap{iteration}'.format(iteration=i), predict=False)
            else:
                logging.info('Iteration ap{iteration} has been imaged already, skipping WSClean.'.format(iteration=i))
            statistic_new = measure_statistic(image_name + '_ap{iteration}-MFS-image.fits'.format(iteration=i))
            statistics.append(statistic_new)
            logging.info('All statistics: ' + ', '.join([str(stats) for stats in statistics]))
            logging.info('Current statistic: ' + str(statistic_current))
            logging.info('New statistic: ' + str(statistic_new))
            if (statistic_new / statistic_current) > 1.02:
                # Image quality improved, predict back.
                logging.info('Image quality improved.')
                stall = 0
            elif statistic_new < statistic_current:
                logging.info('Image quality did not improve, stopping amplitude+phase iterations at this solint.')
                logging.info('Best iteration: {:d}'.format(i-1))
                done = True
            elif (statistic_new / statistic_current) < 1.02:
                logging.info('Image improved less than 2%.')
                stall += 1
                if stall > 3:
                    logging.info('Image improved less than 2% three successive times, selfcal is likely stalling.')
                    done = True
            logging.info('Predicting back to MODEL_DATA.')
            run_wsclean(ms, j=32, mem=75, scale='0.025asec', baseline_averaging=290, no_update_model_required=True, niter=5000, auto_mask=5, local_rms=True, auto_threshold=1.5, weight='briggs', robust=-1, channels_out=12, join_channels=True, stop_negative=True, size=1024, minuvw_m=30000, colname='CORRECTED_DATA_PHASE_ONLY', name=image_name + '_ap{iteration}'.format(iteration=i), predict=True)
            subprocess.call('touch {ms}_has_predict_ampphase_ap{iteration}'.format(ms=ms, iteration=i), shell=True)
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('ms', help='Measurement set to self-cal on.')
    parser.add_argument('--solint-phase', nargs='+', type=int, action='append', help='List of solints to iterate phase-only solves on.')
    parser.add_argument('--nchan-phase', nargs='+', type=int, action='append', help='List of nchans to iterate phase-only solves on.')
    parser.add_argument('--solint-amp', nargs='+', type=int, action='append', help='List of solints to iterate amp+phase solves on.')
    parser.add_argument('--nchan-amp', nargs='+', type=int, action='append', help='List of nchans to iterate amp+phase solves on.')
    parser.add_argument('--image-name', type=str, default='image', help='Image basename.')
    parser.add_argument('--colname', type=str, default='DATA', help='Column name for initial image.')
    arguments = parser.parse_args()
    main(**vars(arguments))

