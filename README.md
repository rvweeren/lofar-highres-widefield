# A LOFAR Long Baseline Widefield Imaging Pipeline
This is an imaging pipeline to produce high quality high resolution widefield images of LOFAR observations including the international stations. It aims to be a pipeline similar as what the ddf-pipeline is for the Dutch LOFAR, and uses DDFacet, DP3 and WSClean.

The end result of the pipeline is currently a 1'' map that is direction-independently calibrated. As input, it requires a dataset, a square DS9 region centered on the pointing center, the corresponding LoTSS reduction and solutions towards an infield long-baseline calibrator. Given this, the pipeline will then execute the following steps:

1. Subtract all sources outside of a given region on the center of the field, using the DD solutions from LoTSS.
2. Produce DP3 parsets to split out potential DDE calibrators.

The pipeline is split up in different steps, to accomodate a Grid workflow. Each can be run individually in the genericpipeline framework as

```
genericpipeline.py <step>.parset
```

Requirements
------------
* DDFacet: https://github.com/saopicc/DDFacet
* ddf-pipeline: https://github.com/mhardcastle/ddf-pipeline
* DP3: https://github.com/lofar-astron/DP3
* LoSoTo: https://github.com/revoltek/losoto
* long baseline pipeline: https://github.com/lmorabit/long_baseline_pipeline
* PyBDSF: https://github.com/lofar-astron/PyBDSF
* WSClean: https://sourceforge.net/p/wsclean/wiki/Home/
