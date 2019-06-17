#!/usr/bin/env bash
python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/runwscleanLB.py --imsize=512 --robust=-1.0 --niter=1000 --longbaseline --pixelscale=0.025 --channelsout=6 -i imagecal --no-tec --uvmin=40000 --lb-nchan-phase=5 --lb-nchan-ap=20 --lb-solint-ap=80 --lb-solint-phase=2 --skymodel=SL333880.skymodel SL333880.ms
