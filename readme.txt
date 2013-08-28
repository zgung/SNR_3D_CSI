The function SNR_with_images is reading SIEMENS 3D CSI dicom and 3D dicom 
images measured in transversal plane.

 !!!!!!!!!!!!!!!!!! readme !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 for working you need my other function called read_ascconv_lenk.m 
       for reading parameters from the dicom
       and a function called bf.m for baseline correction
       http://www.mathworks.com/matlabcentral/fileexchange/24916-baseline-fit

 directory = '~/Patient_name/' - where is a directory called "Spec" with 
       a dicom 3D CSI file and a directory called "Dixon_1.0iso_PAT2_v2_W"
       for water DIXON images

 cho_ppm = 3.2 - exact position of Choline peak if the spectrum is set
       to begin at 8.76 ppm and ends at 0.64, with 1000 Hz bandwidth

 bdwtd = 50 - bandwidth of the peak

 trnct = number of points to truncate at the end of FID

 control = 1 - at first you need to control if the baseline correction 
    is working properly and everything is set all right, if you do not
    wish to continue controling and you need to stop the script running,
    try pressing "ctrl + c" that should stop it

 the output is maximal, mean value of all SNRs of Cho and a table 
 with all SNRs in one row with coordinates, all saved in txt files in Spec directory

 the script also saves the Cho ppm and bandwith that was set for the results
