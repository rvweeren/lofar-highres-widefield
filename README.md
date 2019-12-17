# A LOFAR Long Baseline Widefield Imaging Pipeline

**The pipeline is currently being rewritten in a generic pipeline format and a more parallel fashion oriented on running on clusters this master branch is not maintained at the moment.**

This is an imaging pipeline to produce high quality high resolution widefield images of LOFAR observations including the international stations. It aims to be a pipeline similar as what the ddf-pipeline is for the Dutch LOFAR, and uses DDFacet, DP3 and WSClean.

The end result of the pipeline is currently a 1'' map that is direction-independently calibrated. As input, it requires a dataset, a square DS9 region centered on the pointing center, the corresponding LoTSS reduction and solutions towards an infield long-baseline calibrator. Given this, the pipeline will then:

1. Apply the LoTSS solutions to CS/RS, to arrive at the correct DI calibrated data.
2. Subtract all sources outside of a given region on the center of the field, using the DD solutions from LoTSS.
3. Apply the infield calibrator solutions to DI-correct the international stations.
4. Image the central region of the field at 1'' angular resolution.
5. Produce a PyBDSF catalogue of the image.
6. Produce DP3 parsets to split out potential DDE calibrators.

To run the pipeline, set the appropriate settings in the config file and run:
```
wifi.py lb_widefield.cfg
```

Tyipcal runtime on a Leiden node is about a week for the 1'' DI image.

Requirements
------------
* DDFacet: https://github.com/saopicc/DDFacet
* ddf-pipeline: https://github.com/mhardcastle/ddf-pipeline
* DP3: https://github.com/lofar-astron/DP3
* LoSoTo: https://github.com/revoltek/losoto
* long baseline pipeline: https://github.com/lmorabit/long_baseline_pipeline
* PyBDSF: https://github.com/lofar-astron/PyBDSF
* WSClean: https://sourceforge.net/p/wsclean/wiki/Home/
