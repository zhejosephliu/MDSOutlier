/*
 * mds.cpp
 *
 *  Created on: Jul 24, 2010
 *      Author: Zhe Liu
 */

#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>

#include "head.h"

#include <gsl/gsl_linalg.h>

using namespace std;

void MDS(Mat & ibs, Mat & coordinate, int numofcomp) {

	int i, j;

	if(ibs.size() != ibs[0].size()) {
		cout<<"ERROR!"<<endl;
		cout<<"Numbers of row and column are not equal."<<endl;
		exit(1);
	}

	int length = ibs.size();

	Mat b(length, Vec(length, 0));

	Vec rowmean(length, 0);
	Vec colmean(length, 0);

	double mean = 0;

	for(i = 0; i < length; ++i) {
		for(j = 0; j < length; ++j) {
			rowmean[i] += pow(ibs[i][j], 2);
			colmean[j] += pow(ibs[j][i], 2);
		}
		mean += rowmean[i];
	}

	for(i = 0; i < length; ++i)
		rowmean[i] = rowmean[i] / length;
	for(j = 0; j < length; ++j)
		colmean[j] = colmean[j] / length;

	mean = mean / (length * length);

	for(i = 0; i < length; ++i)
		for(j = 0; j < length; ++j)
			b[i][j] = (-0.5) * (pow(ibs[i][j], 2) - rowmean[i] - colmean[j] + mean);

	// convert into gsl_matrix format
	gsl_matrix * m = gsl_matrix_alloc (length, length);
	gsl_matrix * v = gsl_matrix_alloc (length, length);
	gsl_vector * s = gsl_vector_alloc (length);
	gsl_vector * work = gsl_vector_alloc (length);

	// initialization
	for(i = 0; i < length; ++i) {
		gsl_vector_set(s, i, 0);
		gsl_vector_set(work, i, 0);
		for(j = 0; j < length; ++j) {
			gsl_matrix_set(m, i, j, b[i][j]);
			gsl_matrix_set(v, i, j, 0);
		}
	}

	// call gls_SVD program
	gsl_linalg_SV_decomp(m, v, s, work);

	for(j = 0; j < numofcomp; ++j)
		for(i = 0; i < length; ++i)
			coordinate[i][j] = gsl_matrix_get(m, i, j) * sqrt(gsl_vector_get(s, j));

} ///:~
