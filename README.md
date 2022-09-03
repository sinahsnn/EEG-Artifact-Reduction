# EEG-Artifact-Reduction
In this repository, some blind and semi-blind source separation methods, such as PCA, GEVD, ICA, and DSS, have been implemented to eliminate the EEG artifacts.
## Question 1  
- Ex1.mat data was scattered.
- the PCA algorithm was implemented on the data, and the principal components were extracted. Moreover, the whitened data was plotted. 
- PCA algorithm was implemented with the Matlab function, and the result was similar to the results obtained in the last part.
- SVD decomposition was used, and checked the connections of PCA and SVD
## Question 2  
- Ex2.m contains non-epileptic simulated data, which we contaminated with background EEG signal and muscle noise.
- Sources were extracted using PCA and Com2 (one of the ICA algorithms), and undesirable sources were eliminated.
- Denoised data was transferred to sensor space. Furthermore, RRMSE was calculated for the two methods with different SNRs. 
## Question 3
- Ex3.m includes seizure/non-seizure periods of recorded EEG signals from epileptic patients.
- Main artifacts were detected, and Com2 was applied to detect sources.
- Frequency and spatial characteristics (topo plots) were obtained.
- Undesired sources were removed based on the temporal, spatial, and frequency characteristics, and the denoised signal was constructed.
## Question 4  
the data contains eight channel signals whose duration is 100 seconds. The data is a linear combination of some sources and noises. Three sources are signals with the following characteristics:
* s1(t) : s1 is a periodic triangular signal. 
* s2(t): s2 is a non-stationary signal which was on in some particular intervals.
* s3(t): s3 is a Band limited signal. 
* x(t) = x1(t) + x2(t) + x3(t) + x4(t) :x1 and x2 and x3 are related to s1 , s2 , and s3' s effect respectively and x4 is related to noise sources effect. 
in this part, we used two semi-blind source separation methods named GEVD and DSS for denoising the EEG signals:
- in the first section GEVD method was applied for Periodicity Maximization to extract s1, which is a periodic triangular signal. Moreover, the DSS method was applied. The function f, which is used for denoising in the expectation step, was considered to epoching the data with 400 samples, calculate the mean of data and concatenate the calculated mean with itself to generate data whose samples are similar to the main data. (knowing the period)
- not knowing the exact period time ( we know it is between 3-7 seconds), we applied the described methods in the previous part to denoise the signals. 
- in the first section GEVD method was applied for Nonstationary Maximization to extract s2. Moreover, the DSS method was applied. The function f, which is used for denoising in the expectation step, was considered to make the data zero in the intervals that we know the data is zero. 
- not knowing all the intervals that the signal is on, we applied the described methods in the previous part to denoise the signals. 
- in the first section GEVD method was applied for Spectral Contrast Maximization to extract s3. Moreover, the DSS method was applied. The function f, which is used for denoising in the expectation step, was considered to be a bandpass filter.[(10-15 )Hz]
- not knowing all the exact frequency intervals of signal, we applied the described methods in the last part to denoise the signals. 
## Question 5  
- pure.mat contains a 19-channel EEG signal of a subject with closed eyes. 
- contaminated.m includes an artificially contaminated version of pure. mat. 
- The primary artifact is an electrooculogram signal that is mixed up with an EEG signal. 
- Intervals affected by EOG were detected. 
- Sources were extracted using GEVD and DSS. Undesirable sources were eliminated, and a clean signal was constructed. 




