# Final one to store imaging weights after deciding best resolution.
time wsclean -j 48 -mem 90 -data-column CORRECTED_DATA_LB_AMPPHASE -scale 0.20asec -no-reorder -no-update-model-required -niter 0 -weight briggs -1 -size 1024 1024 -minuvw-m 5000 -taper-gaussian 1asec -store-imaging-weights -channels-out 6 -name wsclean_taper1asec_uvcut5km_230SB_1k_ampphase *.msdpppconcat
#time wsclean -j 48 -mem 90 -data-column CORRECTED_DATA_LB_PHASE -scale 0.20asec -no-update-model-required -niter 5000 -weight briggs -1 -size 10000 10000 -minuvw-m 30000 -taper-gaussian 1asec -baseline-averaging 200 -name wsclean_30km_1asec_taper_10SB_10k_phaseselfcal *.ms
#time wsclean -j 48 -mem 75 -data-column CORRECTED_DATA_LB_AMPPHASE -scale 0.20asec -no-update-model-required -niter 5000 -weight briggs -1 -size 10000 10000 -minuvw-m 30000 -taper-gaussian 1asec -baseline-averaging 200 -name wsclean_30km_1asec_taper_10SB_10k_ampphaseselfcal *.ms

