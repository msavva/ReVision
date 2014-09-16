ReVision
========

This repository contains the code used for classification of images in the UIST 2011 paper [ReVision: Automated Classification, Analysis and Redesign of Chart Images](http://vis.stanford.edu/papers/revision).

The classification code is implemented in MATLAB and is found in the `code/visClass` directory.  There is an external dependency on the vlfeat library included in `code/vlfeat`.  Disclaimer: this is research code that is definitely in need of clean up and more documentation. However, it should be functional, and hopefully easy enough to understand.

The main entry point is in [getTxtAll.m](code/visClass/getTxtAll.m) and global settings are in [setOpt.m](code/visClass/setOpt.m).  There you'll find constant initializations and paths to data.  Once you run the entry script, it will loop through all input image data, extract a set of centroid patches of predefined sizes ("rfs") and then compute feature vectors for the images.  The feature vectors are then saved in Weka ARFF format.  We used Weka's SVM implementation (http://www.cs.waikato.ac.nz/ml/weka/) on this data to train our image classification model and run all our classification experiments.
