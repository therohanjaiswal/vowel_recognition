// vowel_recognition.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>	
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#define PI 3.14
using namespace std;

/* This function computes DC shift*/
float get_dc_shift() 
{
	FILE* file = fopen("dc_shift.txt", "r");	// opening dc_shift.txt in reading mode
	char chars[100] = {'\0'};					// read sample values in characters from file			
	int data;									// store the char sample value in int 
	int sample_count = 0;						// store the total number of sample values
	float sum = 0.0;							// store summation of all data values 
		
	while(!feof(file)) 
	{
		fgets(chars, 100, file);			// reading current line 
		if(!isalpha(chars[0]))				// if the oth character of first line is not an alphabet
		{	
			data = atoi(chars);					
			sum += data;					// summing up sample values
			++sample_count;					// incrementing sample_count						
		}
	}
	fclose(file);			// closing the file
	return sum/sample_count;	// returning average dc shift value
}

/* This function normalize all the sample values and store it into sample array*/
void compute_normalized_data(FILE *file, double *sample_array, float norm_factor, float dc_shift) 
{
	char chars[100] = {'\0'};	// read values in characters from file
	double data = 0;			// data store the char value in int
	int i = 0;					// index to iterate sample_array
	while(!feof(file))			// iterate upto end of file
	{
		fgets(chars, 100, file);
		if(!isalpha(chars[0])) 
		{
			data = atof(chars) - dc_shift;		// subtract dc shift from each sample value
			data *= norm_factor;				// normalize each sample value
			*(sample_array + i) = data;			// store normalized data into sample_array
			i++;								// increment index
		}
	}
}

/* This function computes returns the index of maximum of all the normalized sample values
	present in sample_data array */
double get_max_normalized_data_index(double *sample_data, int sample_count) 
{
	int i = 0;			// iterator
	int result = 0;		// store the index of maximum of all the normalized sample values
	double max = 0;		// store maximum of all the normalized data
	for(i = 0; i < sample_count; i++) 
	{
		if(*(sample_data + i) > max) 
		{
			max = *(sample_data + i);
			result = i;
		}
	}
	return result;		// returns result
}

/* This function computes the Ci values for each frame of training files by using 
   standard formula and store it into Ci[][][] array */
void compute_Ci_for_each_frame(double Ai[13][13], double Ci[10][5][13], int p, int vowel_index, int frame_number)
{
	int m = 1, k = 1, q = 12;	
	double sum = 0;
	for(m = 1; m<=p; m++) {
		sum = 0;
		for(k = 1; k<= m-1; k++) {
			sum += (k * Ci[vowel_index][frame_number][k] * Ai[p][m-k]);
		}
		Ci[vowel_index][frame_number][m] = (Ai[p][m] + sum/m);
	}

	// applying sine window
	 for(m = 1; m<=p; m++) {
		 Ci[vowel_index][frame_number][m] *=  1+((q/2)*sin((3.14*m)/q));
	 }
}

/* This function computes the Ai values for each frame by using standard formula and store it into Ai array*/
void compute_Ai_for_each_frame(double *Ri, double Ai[13][13], double *Ki, double *Ei, int p) 
{
	int i = 1, j = 1;
	double sum = 0;
	*Ei = *Ri;
	for(i = 1; i<= p; i++) {
		sum = 0;
		for(j = 1; j<= i-1; j++) {
			sum += Ai[i-1][j] * *(Ri + i - j);
		}

		*(Ki + i) = (*(Ri + i) - sum)/ *(Ei + i - 1);
		Ai[i][i] = *(Ki + i);

		for(j = 1; j<=i-1; j++) {
			Ai[i][j] = Ai[i-1][j] - *(Ki + i) * Ai[i-1][i-j];
		}

		*(Ei + i) = (1 - *(Ki + i) * *(Ki + i)) * *(Ei + i - 1);
	}
}

/* This function computes the Ri values for each frame by using standard formula and store it into Ri array*/
void compute_Ri_for_each_frame(double *data_array, double *Ri, int p, int sample_size) 
{
	int k = 0, i = 0;
	for(k = 0; k<=p; k++) 
	{
		Ri[k] = 0.0;
		for(i = 1; i<=sample_size-k; i++) 
		{
			*(Ri+k) += *(data_array + i) * *(data_array + i + k);
		}
	}	
}

/* This function computes Ri, Ai and Ci values of each frames for each training file */
void compute_Ri_Ai_Ci_for_training_frame(double Ci[10][5][13], int vowel_index, int frame_index) 
{
	int sample_size = 320;	// frame size
	int p = 12;				// p = 12
	double *data_array = new double[sample_size+1];	// data array to store the sample values of current frame i.e., frame[frame_index]
	double *Ri = new double[p+1];					// store Ri values of current frame

	// For Ai_values	
	double *Ki = new double[p+1];					// store Ki values
	double *Ei = new double[p+1];					// store Ei values
	double Ai[13][13];								// store Ai values	

	ostringstream oss;
	oss << "frame" << frame_index << ".txt";
	string st = oss.str();
    const char *frame_name = st.c_str();
	FILE* file = fopen(frame_name, "r");			// open current frame file in reading mode
	char mystring[150];								// to read lines from file
	double data = 0;								// to store sample values
	int i = 1;										// iterator to store 320 sample values into data array

	// writing on data array
	for(i = 1; i<=sample_size;i++) 
	{
		fgets(mystring, 150, file);					// reading line from current frame file
		data = atof(mystring);						// converting it into float
		*(data_array + i) = data;					// storing into data array
	}
	fclose(file);									// closing current frame file

	compute_Ri_for_each_frame(data_array, Ri, p, sample_size);		// computing Ri values for current frame
	compute_Ai_for_each_frame(Ri, Ai, Ki, Ei, p);					// computing Ai values for current frame
	compute_Ci_for_each_frame(Ai, Ci, p, vowel_index, frame_index);	//computing Ci values for current frame
}

/* This function create frame of each each training file */
void make_frames(double *sample_array, int start_marker) 
{
	int frame_index = 0;			// represents frame number 
	int start = start_marker;		// start index of current frame
	int end = start_marker + 319;	// end index of current frame

	// create five frames (frame0 to frame4)
	for(frame_index = 0; frame_index < 5; frame_index++) 
	{
		start = start_marker;		// start index of current frame
		end = start_marker + 319;	// end index of current frame
		ostringstream oss;
		oss << "frame" << frame_index << ".txt";
		string st = oss.str();
		const char *frame_name = st.c_str();
		FILE *file = fopen(frame_name, "w");			// creating current frame file
		for(int i = start; i<= end; i++) 
			fprintf(file,"%lf\n", *(sample_array + i));	// writing sample values into current frame file 
		fclose(file);
		start_marker = end + 1;		// shifting start_marker to next frame
	}
}

/* This function computes the Ci values for each frame of testing files by using 
   standard formula and store it into Ci[][] array*/
void find_Ci_for_each_test_frame(double Ai[13][13], double Ci[5][13], int p, int frame_number) {
	int m = 1, k = 1, q = 12; 
	double sum = 0;
	for(m = 1; m<=p; m++) {
		sum = 0;
		for(k = 1; k<= m-1; k++) {
			sum += (k * Ci[frame_number][k] * Ai[p][m-k]);
		}
		Ci[frame_number][m] = Ai[p][m] + sum/m;
	}

	// sine window
	for(m = 1; m<=p; m++) {
		 Ci[frame_number][m] *=  1+((q/2)*sin((3.14*m)/q));
	 }
}


/* This function computes Ri, Ai and Ci values of each frames for each testing file */
void compute_Ri_Ai_Ci_for_test_frame(double Ci[5][13], int frame_number) 
{	
	int sample_size = 320;				// denotes frame size of 320 sample values
	int p = 12;							// p = 12
	double *data_array = new double[sample_size+1];	// data array to store the sample values of current frame i.e., frame[frame_index]
	double *Ri = new double[p+1];					// store Ri values of current frame

	// For Ai_values	
	double *Ki = new double[p+1];					// store Ki values
	double *Ei = new double[p+1];					// store Ei values
	double Ai[13][13];								// store Ai values	

	ostringstream oss;
	oss << "frame" << frame_number << ".txt";
	string st = oss.str();
    const char *frame_name = st.c_str();
	FILE* file = fopen(frame_name, "r");			// open current frame file in reading mode
	char mystring[150];								// to read lines from file
	double data = 0;								// to store sample values
	int i = 1;										// iterator to store 320 sample values into data array

	// writing on data array
	for(i = 1; i <= sample_size; i++) {
		fgets(mystring, 150, file);					// reading line from current frame file
		data = atof(mystring);						// converting it into float
		*(data_array + i) = data;					// storing into data array                 
	}
	fclose(file);									// closing current frame file

	compute_Ri_for_each_frame(data_array, Ri, p, sample_size);		// computing Ri values for current frame
	compute_Ai_for_each_frame(Ri, Ai, Ki, Ei, p);					// computing Ai values for current frame
	find_Ci_for_each_test_frame(Ai, Ci, p, frame_number);			// computing Ci values for current frame
}

/* This function finds tokhura distance between test file and one vowel's reference file */
double find_tokhuras_distance(FILE *test_fp, FILE *ref_fp) {
    double weight[12] = { 1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};	// Tokhhura weights 
    int i = 0;							// iterator
	int frame_index = 0;				// represents frame number
    double tokhura_distance = 0;		// store Tokhura distance considering 12 c(r) values of one frame at a time 
	float avg_tokhura_distance = 0;		// store avg tokhura distance computed by taking tokhura distance each frame and divided by 5
    char chars1[150] = {'\0'};			// to store each line from test file
    char chars2[150] = {'\0'};			// to store each line from reference file
    double ct = 0;						// c(t) value
    double cr = 0;						// c(r) value

	// looping through all the 5 frames
    for(frame_index = 0; frame_index < 5; frame_index++) {
        tokhura_distance = 0;			// Tokhura distance of current frame considering 12 c(r) values 
		// reading 12 c(t) and 12 c(r) values from current frame 
		for(i = 0; i < 12; i++) {		
            fgets(chars1, 150, test_fp);					// read line from test file
            ct = atof(chars1);								// assign c(t) value
            fgets(chars2, 150, ref_fp);						// read line from reference file
            cr = atof(chars2);								// assign c(r) value
            tokhura_distance += weight[i]*pow(ct - cr, 2);	// compute Tokhura Distance by standard formula
        }
		// compute average Tokhura Distance by adding Tokhura distance(consdering 12 c values) of each frame and diving by 5
        avg_tokhura_distance += tokhura_distance;			 
    }
    return avg_tokhura_distance/5;	// return average tokhura distance
}

/* This function identify the vowel of testing files by taking minimum of Tokhura Distance
	calculated from all 5 vowel's reference files formed during training */
string identify_vowel() {
    string vowels[5] = {"a", "e", "i", "o", "u"};	// store vowels
	ostringstream oss;
    oss << "test.txt";
    string str = oss.str();
    const char *test_file = str.c_str();
    FILE *test_vowel_fp = fopen(test_file, "r");	// open the test file in reading mode
    
	double min_distance = INT_MAX;					// store the min of tokhura distance of all 5 vowel's reference files formed during training
	double tokhura_dist = 0.0;						// stores tokhura distance	
    string predicted_vowel = "";					// store the predicted vowel
    int i = 0;										// iterator 

	// iterate through all 5 vowel's reference files formed during training
    for(i = 0; i<5; i++) {
        oss.str("");
        oss.clear();
        oss << "ref_" << vowels[i] << ".txt";		 
        str = oss.str();
        const char *ref_file = str.c_str();			
        FILE *ref_vowel_fp = fopen(ref_file,"r");		// open ith vowel's reference file

        tokhura_dist = find_tokhuras_distance(test_vowel_fp, ref_vowel_fp);		// find tokhura distance between test file and ith vowel's reference file 
        printf("Tokhura distance of %s is %lf\n", vowels[i].c_str(), tokhura_dist);	// print tokhura distance calculated
        
		if(tokhura_dist < min_distance) {
            min_distance = tokhura_dist;	// update the min_distance with vowel having minimum Tokhura distance 
            predicted_vowel = vowels[i];	// update the predicted_vowel with vowel having minimum Tokhura distance 
        }
        fclose(ref_vowel_fp);	// close ref_vowel_fp
        rewind(test_vowel_fp);	// rewind test_vowel_fp
    }
    fclose(test_vowel_fp);		// close ref_vowel_fp
    printf("\nPredicted Vowel is - %s \n\n", predicted_vowel.c_str());	// print the predicted vowel
	return predicted_vowel.c_str();		// return the predicted vowel
}	

/* This function performs the testing of all the 10 testing files of each vowel*/
void testing_of_one_vowel(int vowel_index) {
	string vowels[5] = {"a", "e", "i", "o", "u"};	// store vowels
	int i = 0, j = 0;				// iterators
	int frame_index = 0;			// iterator representing frame number of each frames					
	int correct_predictions = 0;	// represents the no of testing files of current vowel i.e., vowels[vowel_index] is predicted correct
	string predicted_vowel = "";	// represents the predicted vowel from one testing file of current vowel

	// iterate through all the 10 testing files of current vowel i.e., vowels[vowel_index]
	for(i = 11; i <= 20; i++) {
		ostringstream oss;
		oss << "testing_set/214101042_" << vowels[vowel_index] << "_" << i << ".txt";	
		string st = oss.str();
		printf("Testing file taken - %s\n", st.c_str());				// print test file name of current vowel 
		const char *test_file = st.c_str();			
		FILE* test_fp = fopen(test_file, "r");		// open test file in reading mode
		char chars[150];							// read every line from the file

		// reading the ith training file and finding max sample value and total number of samples present 
		int max_sample = 0;					// max_sample stores max of sample value from ith training file
		int my_data = 0;					// to store the lines read from a file
		float norm_factor = 0.0;			// store normalization factor
		float dc_shift = 0.0;				// store dc shift 
		long sample_count = 0;				// store total number of sample values in ith training file
		int start_marker = 0;				// store the index from which reading reading of frames will happen
		int max_normalized_data_index;		// store maximum sample value of normalized data
		double *sample_array;				// store the sample values of ith training file
		double Ci[5][13];					// store all the 12 Ci values for each frame of ith testing file

		dc_shift = get_dc_shift();			// computing dc shift

		// reading the ith training file and finding max sample value and total number of samples present 
		while(!feof(test_fp)) {
			fgets(chars, 100, test_fp);			// reading 1st line from the file
			if(!isalpha(chars[0])) {			// if first read character is not an alphabet
				my_data = atoi(chars);			//convert each char into int					
				max_sample = my_data > max_sample ? my_data : max_sample;	// findind maximum of all sample data
				++sample_count;												// finding total number of sample data
			}
		}

		sample_array = new double[sample_count];	// creating array to store the sample values of ith testing file
		max_sample -= dc_shift;						// eliminating dc shift from max sample value
		norm_factor = 5000.0/max_sample;			// computing normalization factor

		rewind(test_fp);
		compute_normalized_data(test_fp, sample_array, norm_factor, dc_shift);	// normalizing sample values of ith testing file
		fclose(test_fp);

		max_normalized_data_index = get_max_normalized_data_index(sample_array, sample_count);	// getting index of max normalized sample value 
		start_marker = max_normalized_data_index - 800;	// frames will be considered from the 800 samples before the max normalized data index

		// make 5 frames(frame0 to frame4) each of size 320 for ith testing file
		make_frames(sample_array, start_marker);
	
		// compute Ri, Ai and Ci values for each frames for ith testing file
		for(frame_index = 0; frame_index < 5; frame_index++){
			compute_Ri_Ai_Ci_for_test_frame(Ci, frame_index);
		}

		oss.str("");
		oss.clear();
		oss << "test.txt";
		st = oss.str();
		test_file = st.c_str();		
		test_fp = fopen(test_file, "w");		// creating test file for ith testing file to store Ci values of it 5 frames(frame0 to frame4)
		
		// looping through all the frames
		for(frame_index = 0; frame_index < 5; frame_index++) {
			for(j = 1; j<13; j++) {
				fprintf(test_fp, "%lf\n", Ci[frame_index][j]);	// storing all the 12 Ci values for each frame in Ci array
			}
		}
		fclose(test_fp);						// closing file pointer
		predicted_vowel = identify_vowel();		// identify the vowel of ith testing file of current vowel i.e., vowels[vowel_index]
		if(predicted_vowel.compare(vowels[vowel_index].c_str()) == 0)
			++correct_predictions;
	}
	printf("Accuracy : %d percent\n", correct_predictions * 10);
}

/* This function performs the training of all the 10 training files of each vowel*/
void training_of_one_vowel(int vowel_index) 
{	
	string vowels[] = {"a", "e", "i", "o", "u"};	// store vowels
	ostringstream oss1;
	oss1 << "ref_" << vowels[vowel_index] << ".txt";
	string str = oss1.str();
	const char *ref_file = str.c_str(); 
	FILE *ref_fp = fopen(ref_file,"w");				// file pointer to create reference file of current vowel i.e., vowels[vowel_index]

	double Ci[10][5][13];							// Ci stores the Ci values of all 5 frames of 10 training files of current vowel i.e.,vowels[vowel_index]
	double avg_ci[5][13];							// avg_ci store the average ci values(from c1 to c12) for 5 frames i.e., it contains 60 values
	int i = 0, j = 0, k = 0;						// iterators
	int frame_index = 0;							// iterator representing frame number of each frames						

	// iterate through all the 10 training files of current vowel i.e., vowels[vowel_index]
	for(i = 0; i < 10; i++) {
		ostringstream oss;
		oss << "training_set/214101042_" << vowels[vowel_index] << "_" << i+1 << ".txt";
		string st = oss.str();
		const char *vowel_file = st.c_str();
		FILE* fp = fopen(vowel_file, "r");	// opening ith training file of current vowel in reading mode
		char chars[150];					

		int max_sample = 0; 				// max_sample stores max of sample value from ith training file
		int my_data = 0;					// to store the lines read from a file
		float norm_factor = 0.0;			// store normalization factor
		float dc_shift = 0.0;				// store dc shift 
		long sample_count = 0;				// store total number of sample values in ith training file
		int start_marker = 0;				// store the index from which reading reading of frames will happen
		int max_normalized_data_index = 0;	// store maximum sample value of normalized data
		double *sample_array;				// store the sample values of ith training file
			
		dc_shift = get_dc_shift();			// computing dc shift

		// reading the ith training file and finding max sample value and total number of samples present 
		while(!feof(fp)) {
			fgets(chars, 100, fp);					// reading 1st line from the file
			if(!isalpha(chars[0])) {				// if first read character is not an alphabet
				my_data = atoi(chars);				//convert each char into int
				max_sample = my_data > max_sample ? my_data : max_sample;	// findind maximum of all sample data
				++sample_count;						// finding total number of sample data
			}
		}

		sample_array = new double[sample_count];	// creating array to store the sample values of ith training file
		max_sample -= dc_shift;						// eliminating dc shift from max sample value
		norm_factor = 5000.0/max_sample;			// computing normalization factor

		rewind(fp);
		compute_normalized_data(fp, sample_array, norm_factor, dc_shift);	// normalizing sample values of ith training file
		fclose(fp);

		max_normalized_data_index = get_max_normalized_data_index(sample_array, sample_count);	// getting index of max normalized sample value 
		start_marker = max_normalized_data_index - 800;	// frames will be considered from the 800 samples before the max normalized data index

		// make 5 frames(frame0 to frame4) each of size 320 for ith training file
		make_frames(sample_array, start_marker);

		// compute Ri, Ai and Ci values for each frames for ith training file
		for(frame_index = 0; frame_index < 5; frame_index++){
			compute_Ri_Ai_Ci_for_training_frame(Ci, i, frame_index);
		}
	}

	// making one reference file containing avg ci values from Ci array for current vowel
	double sum = 0;
	for(i = 0; i < 5; i++) {
		for(j = 1; j < 13; j++) {
			sum = 0;
			for(k = 0; k < 10; k++) {
				sum += Ci[k][i][j];
			}
			avg_ci[i][j] = sum/10;
			fprintf(ref_fp, "%lf\n", avg_ci[i][j]);
		}
	}
	fclose(ref_fp);
}

int _tmain(int argc, _TCHAR* argv[])
{	
	// ***************** Training of each vowel ***************************					
	string vowels[5] = {"a", "e", "i", "o", "u"};	// store vowels
	int vowel_index = 0;	// iterator representing index of current vowel

	// looping for training of each vowel i.e., making reference files of each vowel 
	for (vowel_index = 0; vowel_index < 5; vowel_index++) {
		training_of_one_vowel(vowel_index);	// training of ith vowel from it's 10 training files	
	}

	// ***************** Testing of each vowel ***************************
	for(vowel_index = 0; vowel_index < 5; vowel_index++) {
		printf("\n***************** Testing of vowel %s ***************************\n\n", vowels[vowel_index].c_str());
		testing_of_one_vowel(vowel_index);
	}

	printf("Program Termination!");
	return 0;
}


