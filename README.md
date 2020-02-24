# A LOFAR Long Baseline Widefield Imaging Pipeline
This is an imaging pipeline to produce high quality high resolution widefield images of LOFAR observations including the international stations. It aims to be a pipeline similar as what the ddf-pipeline is for the Dutch LOFAR, and uses DDFacet, DP3 and WSClean.

The end result of the pipeline is currently a 1'' map that is direction-independently calibrated. As input, it requires a dataset, a square DS9 region centered on the pointing center, the corresponding LoTSS reduction and solutions towards an infield long-baseline calibrator. Given this, the pipeline will then execute the following steps:

1. Subtract all sources outside of a given region on the center of the field, using the DD solutions from LoTSS.
2. Produce DP3 parsets to split out potential DDE calibrators.

The pipeline is split up in different steps, to accomodate a Grid workflow. Each can be run individually in the genericpipeline framework as

```
genericpipeline.py <step>.parset
```

It is mainly designed to work under the GRID_LRT package. This means parallelization is handled on a higher level, i.e. when a job runs on individual subbands, as many jobs as subbands will/should be submitted, running the parsets on that.

Durations are measured on the Spider cluster. This has approximately 15 nodes of 32 cores and NVMe SSDs capable of 6 GB/s read/write as local scratch.

Steps in the pipeline are (or will be):

1. Subtract the 6" LoTSS map from the input data. This can operate on the whole bandwidth at once. Duration: <~24 hours
2. Find DDE calibrator candidates from the LoTSS catalog and split them out. These are sources brighter than 10 mJy/beam peak flux. This operates on the 10SB blocks individually. Duration: ~29 hours on NVMe SSDs, 24 simultaneous 8 core jobs.
3. Collect all SB per source and selfcal on the DDE calibrator candidates. This operates on the full bandwidth of each calibrator individually. Duration: 6 hours, 162 jobs.
-- untested below this --
4. Image the field at approximately 1.2" resolution.
5. Divide the field in 9x9 facets and created 1" map-subtracted datasets for each of them.
6. Image each of the facets.

Requirements
------------
* DDFacet: https://github.com/saopicc/DDFacet
* ddf-pipeline: https://github.com/mhardcastle/ddf-pipeline
* DP3: https://github.com/lofar-astron/DP3
* LoSoTo: https://github.com/revoltek/losoto
* long baseline pipeline: https://github.com/lmorabit/long_baseline_pipeline
* PyBDSF: https://github.com/lofar-astron/PyBDSF
* WSClean: https://sourceforge.net/p/wsclean/wiki/Home/
