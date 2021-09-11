// vowel_recognition.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>	
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

void find_Ri_for_each_frame(double *data_array, double *Ri, int p, int sample_size) {
	for(int k = 0; k<=p; k++) {
		Ri[k] = 0.0;
		for(int i = 1; i<=sample_size-k; i++) {
			*(Ri+k) += *(data_array + i) * *(data_array + i + k);
		}
	}	
}

void find_Ai_for_each_frame(double *Ri, double Ai[13][13], double *Ki, double *Ei, int p) {
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

void find_Ci_for_each_frame(double Ai[13][13], double Ci[10][5][13], int p, int vowel_index, int frame_number) {

	int m = 1, k = 1; 
	double sum = 0;
	for(m = 1; m<=p; m++) {
		sum = 0;
		for(k = 1; k<= m-1; k++) {
			sum += (k * Ci[vowel_index][frame_number][k] * Ai[p][m-k]);
		}
		Ci[vowel_index][frame_number][m] = Ai[p][m] + sum/m;
	}

	// for(int i = 1; i<=p; i++) {
	// 	printf("%f\n", Ci[vowel_index][frame_number][i]);
	// }
	// printf("---------------------------------\n");
}

float get_dc_shift() 
{
	FILE* file = fopen("dc_shift.txt", "r");	// opening dc_shift.txt in r mode
	char chars[100] = {'\0'};					// read sample values in characters from file			
	int data, sampleCount = 0;					// data store the char sample value in int, sampleCount store the total number of samples
	float sum = 0.0;							// sum stores summation of all data values 
		
	while(!feof(file)) 
	{
		fgets(chars, 100, file);
		if(!isalpha(chars[0])) 
		{
			data = atoi(chars);
			sum += data;
			++sampleCount;
		}
	}
	fclose(file);			// closing the file
	return sum/sampleCount;	// returning average dc shift value
}

void compute_normalized_data(FILE *file, double *sample_array, float norm_factor, float dc_shift) 
{
	char chars[100] = {'\0'};	// read values in characters from file
	double data = 0;				// data store the char value in int
	int i = 0;

	while(!feof(file)) 
	{
		fgets(chars, 100, file);
		if(!isalpha(chars[0])) 
		{
			data = atof(chars) - dc_shift;		// subtract dc shift from each sample value
			data *= norm_factor;					// normalize each sample value
			*(sample_array + i) = data;
			// fprintf(outputFile,"%lf -- %d\n",data, i);	// writing normalized data into output.txt
			i++;
		}
	}
	// fclose(outputFile);	// closing the file
}

double get_max_normalized_data_index(double *sample_data, int sample_count) {
	int i = 0, index = 0;
	double max = 0;
	for(i = 0; i < sample_count; i++) {
		if(*(sample_data + i) > max) {
			max = *(sample_data + i);
			index = i;
		}
	}
	return index;
}

void make_frames(double *sample_array, int start_marker, int frame_number) {
	int start = start_marker;
	int end = start_marker + 319;
	// printf("start marker - %d , end marker - %d\n", start, end);
	ostringstream oss;
	oss << "frame" << frame_number << ".txt";
	string st = oss.str();
    const char *frame_name = st.c_str();
	FILE *file = fopen(frame_name, "w");
	for(int i = start; i<= end; i++) {
		fprintf(file,"%lf\n", *(sample_array + i));
	}
	fclose(file);
}

void compute_Ci(double Ci[10][5][13], int vowel_index, int frame_number) {
	
	int sample_size = 320;
	int p = 12;
	double *data_array = new double[sample_size+1];
	double *Ri = new double[p+1];

	// For Ai_values
	double *Ki = new double[p+1];
	double *Ei = new double[p+1];
	double Ai[13][13];

	ostringstream oss;
	oss << "frame" << frame_number << ".txt";
	string st = oss.str();
    const char *frame_name = st.c_str();
	FILE* file = fopen(frame_name, "r");
	char mystring[150];
	double data = 0;
	int i = 1;

	for(i = 1; i<=sample_size;i++) {
		fgets(mystring, 150, file);
		data = atof(mystring);
		*(data_array + i) = data;
	}
	fclose(file);

	find_Ri_for_each_frame(data_array, Ri, p, sample_size);
	find_Ai_for_each_frame(Ri, Ai, Ki, Ei, p);
	find_Ci_for_each_frame(Ai, Ci, p, vowel_index, frame_number);
}

void find_Ci_for_each_test_frame(double Ai[13][13], double Ci[5][13], int p, int frame_number) {
	int m = 1, k = 1; 
	double sum = 0;
	for(m = 1; m<=p; m++) {
		sum = 0;
		for(k = 1; k<= m-1; k++) {
			sum += (k * Ci[frame_number][k] * Ai[p][m-k]);
		}
		Ci[frame_number][m] = Ai[p][m] + sum/m;
	}
}

void compute_Ai_Ri_Ci_for_test_frame(double Ci[5][13], int frame_number) {
	
	int sample_size = 320;
	int p = 12;
	double *data_array = new double[sample_size+1];
	double *Ri = new double[p+1];

	// For Ai_values
	double *Ki = new double[p+1];
	double *Ei = new double[p+1];
	double Ai[13][13];

	ostringstream oss;
	oss << "frame" << frame_number << ".txt";
	string st = oss.str();
    const char *frame_name = st.c_str();
	FILE* file = fopen(frame_name, "r");
	char mystring[150];
	double data = 0;
	int i = 1;

	for(i = 1; i<=sample_size;i++) {
		fgets(mystring, 150, file);
		data = atof(mystring);
		*(data_array + i) = data;
	}
	fclose(file);

	find_Ri_for_each_frame(data_array, Ri, p, sample_size);
	find_Ai_for_each_frame(Ri, Ai, Ki, Ei, p);
	find_Ci_for_each_test_frame(Ai, Ci, p, frame_number);
}

double find_tokhuras_distance(FILE *test_fp, FILE *ref_fp) {
    double weight[12] = { 1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
    int i = 0, j = 0;
    double tokhura_distance = 0, avg_tokhura_distance = 0;
    char chars1[150] = {'\0'};
    char chars2[150] = {'\0'};
    double ct = 0;
    double cr = 0;
    for(i = 0; i<5; i++) {
        tokhura_distance = 0;
        for(j = 0; j<12; j++) {
            fgets(chars1, 150, test_fp);
            ct = atof(chars1);
            fgets(chars2, 150, ref_fp);
            cr = atof(chars2);
            tokhura_distance += weight[j]*pow(ct - cr, 2);
        }
        avg_tokhura_distance += tokhura_distance;
    }
    return avg_tokhura_distance/5;
}

void identify_vowel() {
    ostringstream oss;
    oss << "testing_set/test.txt";
    string str1 = oss.str();
    const char *test_file = str1.c_str();
    FILE *test_vowel_fp = fopen(test_file, "r");
    double min_distance = INT_MAX, tokhura_dist = 0.0;
    string result = "";
    string vowels[5] = {"a", "e", "i", "o", "u"};
    int i=0; 
    for(i = 0; i<5; i++) {
        oss.str("");
        oss.clear();
        oss << "testing_set/ref_" << vowels[i] << ".txt";
        string str2 = oss.str();
        const char *ref_file = str2.c_str();
        FILE *ref_fp = fopen(ref_file,"r");
        tokhura_dist = find_tokhuras_distance(test_vowel_fp, ref_fp);
        printf("Tokhura distance of %s is %lf\n", vowels[i].c_str(), tokhura_dist);
        if(tokhura_dist < min_distance) {
            min_distance = tokhura_dist;
            result = vowels[i];
        }
        fclose(ref_fp);
        rewind(test_vowel_fp);
    }
    fclose(test_vowel_fp);
    printf("result string - %s ----------------\n\n", result.c_str());
}

int _tmain(int argc, _TCHAR* argv[])
{	
	// Training ------------------------------------------
	int l = 0, i = 0, j = 0, k = 0;
	string vowels[5] = {"a", "e", "i", "o", "u"};

	for (l = 0; l<5; l++) {
		ostringstream oss1;
		oss1 << "testing_set/ref_" << vowels[l] << ".txt";
		string str = oss1.str();
		const char *ref_file = str.c_str();
		FILE *ref_fp = fopen(ref_file,"w");
		i = 0, j = 0, k = 0;
		double Ci[10][5][13];
		double avg_ci[5][13];
		double sum = 0;

		// for 10 samples of one vowel
		for(i = 0; i < 10; i++) {
			ostringstream oss;
			oss << "training_set/214101042_" << vowels[l] << "_" << i+1 << ".txt";
			string st = oss.str();
			const char *vowel_file = st.c_str();
			FILE* fp = fopen(vowel_file, "r");
			char chars[150];
			int max_sample = 0, my_data = 0;
			float norm_factor = 0.0, dc_shift = 0.0;
			long sample_count = 0;
			int start_marker = 0;
			double *sample_array;

			dc_shift = get_dc_shift();

			while(!feof(fp)) {
				fgets(chars, 100, fp);
				if(!isalpha(chars[0])) {
					my_data = atoi(chars);									//convert each char into int
					max_sample = my_data > max_sample ? my_data : max_sample;	// findind maximum of all sample data
					++sample_count;										// finding total number of sample data
				}
			}

			sample_array = new double[sample_count];
			max_sample -= dc_shift;
			norm_factor = 5000.0/max_sample;

			// printf("norm factor - %f\n", norm_factor);

			rewind(fp);
			compute_normalized_data(fp, sample_array, norm_factor, dc_shift);
			fclose(fp);

			int max_normalized_data_index = get_max_normalized_data_index(sample_array, sample_count);
			//printf("max_normalized_ index ------ %d\n", max_normalized_data_index);
	
			start_marker = max_normalized_data_index - 800;
	
			// j is frame number from 0 to 4
			for(j = 0; j < 5; j++) {
				make_frames(sample_array, start_marker, j);
				start_marker += 320;
			}

			for(j = 0; j < 5; j++){
				compute_Ci(Ci, i, j);
			}

			// printf("\n ---------End of file %d -----------\n", i);
		}
	
		// making ref files for each vowel
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

	// Testing ------------------------------
	i = 0, j = 0, k = 0;
	for(i = 0; i<5; i++) {
		for(j = 11; j<=20; j++) {
			ostringstream oss;
			oss << "testing_set/214101042_" << vowels[i] << "_" << j << ".txt";
			string st = oss.str();
			printf("File taken %s\n", st.c_str());
			const char *test_file = st.c_str();
			FILE* test_fp = fopen(test_file, "r");
			char chars[150];
			int max_sample = 0, my_data = 0;
			float norm_factor = 0.0, dc_shift = 0.0;
			long sample_count = 0;
			int start_marker = 0;
			double *sample_array;
			double Ci[5][13];
			dc_shift = get_dc_shift();
			int i = 0, j = 0;

			while(!feof(test_fp)) {
				fgets(chars, 100, test_fp);
				if(!isalpha(chars[0])) {
					my_data = atoi(chars);									//convert each char into int
					max_sample = my_data > max_sample ? my_data : max_sample;	// findind maximum of all sample data
					++sample_count;										// finding total number of sample data
				}
			}

			sample_array = new double[sample_count];
			max_sample -= dc_shift;
			norm_factor = 5000.0/max_sample;

			// printf("norm factor - %f\n", norm_factor);

			rewind(test_fp);
			compute_normalized_data(test_fp, sample_array, norm_factor, dc_shift);
			fclose(test_fp);

			int max_normalized_data_index = get_max_normalized_data_index(sample_array, sample_count);
			//printf("max_normalized_ index ------ %d\n", max_normalized_data_index);
	
			start_marker = max_normalized_data_index - 800;
	
			// frames 0 to 4
			for(j = 0; j < 5; j++) {
				make_frames(sample_array, start_marker, j);
				start_marker += 320;
			}

			for(j = 0; j < 5; j++){
				compute_Ai_Ri_Ci_for_test_frame(Ci, j);
			}

			oss.str("");
			oss.clear();
			oss << "testing_set/test.txt";
			st = oss.str();
			test_file = st.c_str();
			test_fp = fopen(test_file, "w");

			for(i=0; i<5; i++) {
				for(j = 1; j<13; j++) {
					fprintf(test_fp, "%lf\n", Ci[i][j]);
				}
			}
			fclose(test_fp);
			identify_vowel();
		}
		printf("--------------------------------------------------------------------------------------\n");
	}
	printf("Program Termination!");
	return 0;
}





