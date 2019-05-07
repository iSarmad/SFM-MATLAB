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

* Matlab 2017a (Other Version also might work, this is what I used)
lease refer to the insturctions below to produce results for this assignment.

1. I did not follow the provided template.

2. I made a new one and used data structures such as cell etc.

3. Go to main.m and select 'dataSet' variable 1 for Default dataset and 2 for my own dataset

4. Set ‘highRes’ variable in the main.m to 1 to view extremely high-resolution images. 
(Warning: This will take 20 to 30 minutes to run, but you can do this reproduce figure in this report.).  

5. By default, the flag is set to 0 so that the code runs in 2 minutes but the results will
not be as high quality as the report.

6. Please use the matlab figure to visualize the results instead of the mesh lab

7. Figure files have also been saved in the current directory to visualize without running the high 
resolution version
 File Name : 
  (i). two View : "twoview.fig" for your dataset, and "twoview2.fig" for my data set
   (ii). Growing step: "finalGrowingStep1.fig" for your dataset, "finalGrowingStep2.fig" for my data set
    (iii).Also have a look at SFM1.png and SFM2.png for growing step with your dataset

8. You might need to mex the function for 5 point algorithm again on your computer 

9. Finally to run Bundle Adjustment for two view, set 'optim' variable to 1. By deafult it will be 0.

Thank you.

## Step1 : Feature Extraction and Matching Step



## Step2 : Initialization Step ( Two View Step)


## Step 3 : Growing Step



## License

This project is licensed under the MIT License. 
For specific helper function used in this repository please see the license agreement of the Repo linked in Acknowledgement section
## Acknowledgments
My implementation has been inspired from the following sources.

* [Camera Geometry Matlab](https://uk.mathworks.com/matlabcentral/fileexchange/47032-camera-geometry-algorithms?focused=3822640&tab=function) Mostly just looked at this for sanity check.
