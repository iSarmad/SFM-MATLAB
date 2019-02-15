# SFM-MATLAB
### Teaser 
![alt text](https://github.com/iSarmad/SFM-MATLAB/blob/master/SFM2.png)

Note: This repository is under construction. Thank you for your patience.
This repository is the implementation of "Structure From Motion from Multiple Views" in MATLAB. Refer to following for in depth understanding:

* [Photo Tourism](http://phototour.cs.washington.edu/Photo_Tourism.pdf) Photo Tourism: Exploring Photo Collections in 3D

## About the Dataset
The provided dataset of images was taken from a calibrated camera. Therefore the Camera matrix was known to me. For images provided in this repository, you do not need to perform any camera calibration. But if you want to use your own dataset, first you have to perform camera calibration using the following in detail guide: 

* [Camera Calibration Toolbox](http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/own_calib.html) Doing your own camera calibration


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

What things you need to install:

1. Matlab 2017a (Other Version also might work, this is what I used)


Average run time depends on the number of features used for matching. I made a switch where you can select a high resolution version( takes 20 minutes to run). And a low resolution version (takes 3 to 4 minutes to run)
 
## Step1 : Feature Extraction and Matching Step



## Step2 : Initialization Step ( Two View Step)


## Step 3 : Growing Step



## License

This project is licensed under the MIT License. 
For specific helper function used in this repository please see the license agreement of the Repo linked in Acknowledgement section
## Acknowledgments
My implementation has been inspired from the following sources.

* [Camera Geometry Matlab](https://uk.mathworks.com/matlabcentral/fileexchange/47032-camera-geometry-algorithms?focused=3822640&tab=function) Mostly just looked at this for sanity check.
