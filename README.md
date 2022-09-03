# EEG-Artifact-Reduction
In this repository, some blind and semi-blind source separation methods, such as PCA, GEVD, ICA, and DSS, have been implemented to eliminate the EEG artifacts.
## Question 1  
- Ex1.mat data was scattered.
- the PCA algorithm was implemented on the data, and the principal components were extracted. Moreover, the whitened data was plotted. 
- PCA algorithm was implemented with the Matlab function, and the result was similar to the results obtained in the previous part.
- SVD decomposition was used, and checked the connections of PCA and SVD
## Question 2  
- Contains non-epileptic simulated data, which we contaminated with background EEG signal and muscle noise
- Extracted sources using PCA and Com2 (on of ICA algorithms) and eliminated undesirable sources
- Transferred the data to sensor space. Moreover, calculated RRMSE for the two methods. 
## Question 3  
aksdfjlskdfjsdlkfjsdlkf
## Question 4  
the data contains 8 channels signals whose duration are 100 second. the data is a linear combination of some sources and noises. three sources are signals with the following charactristics:
* s1(t) : s1 is a periodic triangular signal. 
* s2(t) : s2 is a non-statinary signal which was on in some special intervals.
* s3(t) : s3 is a Band limited signal. 
* x(t) = x1(t) + x2(t) + x3(t) + x4(t) 
* x1 and x2 and x3 are related to s1 , s2 , and s3' s effect recpectivelt and x4 is related to noise sources effect. 
* in this part we use two semi bline source speration method named GEVD and DSS for denoising the eeg signals. 
-
-
-
## Question 5  
aksdfjlskdfjsdlkfjsdlkf





