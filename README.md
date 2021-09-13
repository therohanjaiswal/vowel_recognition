
	Vowel Recognition
----------------------------------------------------------

I. Files used to execute the program : <br/>
	* Here, folder training_set contains 10 recorded files of each vowels which is used to generate
	  train the model by generating reference file of each vowel. Folder testing_set contains  10 
	  recorded files of each vowels which is used for vowel prediction. One dc_shift.txt file is also
	  to compute the DC shift. File vowel_recognition.cpp contains the actual cpp code.

II. Idea used to consider 5 frames i.e., frame 0 to frame 4 of each file : <br/>
	 * Here, after normalizing the sample values, I find out the index of maximum sample value of 
	  all sample values present in the file. From the index of maximum sample value, I take 800 
	  previous sample values and 800 next sample values. Overall it forms 1600 sample values which
	  forms 5 frames each containing 1600/5 = 320 sample values.

III. Finding Tokhura Distance : <br/>
	 * Since, each file has 5 frames and each frame has 12 ci values. So, overall each file has
	  12 * 5 = 60 "c(t)/c(r)" values and we can find the Tokhura distance containing 12 "c(t)/c(r)" 
	  values only. So, I computed the Tokhura distance of each frame containing 12 "c(t)/c(r)" values
	  and find out the average Tokhura distance by summing up the Tokhura distance of each frame and 
	  dividing it by 5.

IV. Steps involved in the program: <br/>
    ================ Training ================ <br/>
    1. Take the 10 training files of each vowel present in "training_set" folder. <br/>
    2. Now, perform the following operations for each training file. <br/>
       2.1. Eliminate DC Shift from the sample values. <br/>
       2.2. Take out the maximum sample value and total number of samples. <br/>
       2.3. Compute the normalization factor. <br/>
       2.4. Normalize all the sample values within a range of -5000 to +5000 and store it into an array.  <br/>
       2.5. Find the index of maximum normalized sample value. <br/>
       2.6. From this index, take the previous 800 values and next 800 values which forms 1600 sample values
            required to make 5 frames each consisting of 1600/5 = 320 sample values. <br/>
       2.7. Compute Ri, Ai and Ci values of each frames. <br/>
    3. After finding the Ci values, compute the average of ci values from 10 training files of each vowel. <br/>
    4. Create the reference file containing the average of ci values from 10 training files of each vowel. <br/>
     ================ Testing ================ <br/>
    5. Take the 10 testing files of each vowel present in "testing_set" folder. <br/>
    6. Now, perform the following operations for each testing file. <br/>
       6.1. Eliminate DC Shift from the sample values. <br/>
       6.2. Take out the maximum sample value and total number of samples. <br/>
       6.3. Compute the normalization factor. <br/>
       6.4. Normalize all the sample values within a range of -5000 to +5000 and store it into an array. 
       6.5. Find the index of maximum normalized sample value. <br/>
       6.6. From this index, take the previous 800 values and next 800 values which forms 1600 sample values
            required to make 5 frames each consisting of 1600/5 = 320 sample values. <br/>
       6.7. Compute Ri, Ai and Ci values of each frames. <br/>
       6.8. Find out the Tokhura distance by considering the ci values from this testing file and ci values 
          from reference file of all the vowels. <br/>
       6.9. Reference file of the vowel which generates the minimum Tokhura Distance is the predicted vowel. <br/>
    7. In the end, depending upon the number of correct predicted vowels out of 10 testing files of each vowel,
       it prints the accuracy of the model. <br/>
	
Submitted By : Rohan Jaiswal, MTech CSE, IIT Guwahati. Roll no - 214101042.
