# LOFAR Long Baseline Widefield Imaging Pipeline
This is an imaging pipeline to produce high quality 1'' images of LOFAR observations including the international stations. It aims to be a pipeline similar as what the ddf-pipeline is for the Dutch LOFAR, and uses DDFacet, DP3 and WSClean.

The end result of the pipeline is a 1'' map (final size to be determined) that is direction-independently calibrated. As input, it requires a dataset, a square DS9 region centered on the pointing center, the corresponding LoTSS reduction and solutions towards an infield long-baseline calibrator. Given this, the pipeline will then:

1. apply the LoTSS solutions to CS/RS, to arrive at the correct DI calibrated data.
2. subtract all sources outside of a given region on the center of the field, using the DD solutions from LoTSS.
3. apply the infield calibrator solutions to DI-correct the international stations.
4. image the central region of the field at 1'' angular resolution.

To run the pipeline, set the appropriate settings in the config file and run:
```
wifi.py lb_widefield.cfg
```

Requirements
------------
* DDFacet: https://github.com/saopicc/DDFacet
* ddf-pipeline: https://github.com/mhardcastle/ddf-pipeline
* DP3: https://github.com/lofar-astron/DP3
* LoSoTo: https://github.com/revoltek/losoto
* long baseline pipeline: https://github.com/lmorabit/long\_baseline\_pipeline
* PyBDSF: https://github.com/lofar-astron/PyBDSF
* WSClean: https://sourceforge.net/p/wsclean/wiki/Home/
