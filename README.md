# Blue-Noise-and-Nystrom-extension-for-Change-detection
Graph based changed detection approach presented at IGARSS 2021

Here you will find all the code used to generate the change maps in fourteen cases of study.

bs_cd.m is the main code.
cohensKappa.m will generate all the results showed in the paper and also the change map with respec to FA, MA, P, R, K, OE.

Please if you use the datasets and the codes (including the toolboxes) refer them properly.

## Datasets

To get access to the dataset please go to the next [repository](https://github.com/DavidJimenezS/GBF-CD/tree/master/Data) to download them.

## Rquirements

In order to run the code, you will need the following toolboxes:

* Graph analysis [toolbox](http://leogrady.net/software/) for graph cut segmentation 
* Graph signal processing [toolbox](https://epfl-lts2.github.io/gspbox-html) for graph smoothness prior. 
* Blue-Noise sampling [toolbox](https://github.com/jhonygiraldo/Blue-Noise-Sampling-on-Graphs) for graph cut segmentation.
