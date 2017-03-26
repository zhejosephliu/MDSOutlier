/*
 * print.cpp
 *
 *  Created on: Jul 24, 2010
 *      Author: Zhe Liu
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstring>

#include "head.h"

using namespace std;

void Printvector(const Vec & vector, const char * filename) {

	unsigned int i;
	const char * screen = "screen";

	if(strcmp(filename, screen) == 0) {
		for(i = 0; i < vector.size(); ++i) {
			cout<<vector[i]<<endl;
		}
	}
	else {
		ofstream out(filename);
		for(i = 0; i < vector.size(); ++i) {
			out<<vector[i]<<endl;
		}
		out.close();
	}

} ///:~

void Printmatrix(const Mat & matrix, const char * filename) {

	unsigned int i, j;
	const char * screen = "screen";

	if(strcmp(filename,screen) == 0) {
		for(i = 0; i < matrix.size(); ++i) {
			for(j = 0; j < matrix[0].size(); ++j)
				cout<<matrix[i][j]<<" ";
			cout<<endl;
		}
	}
	else {
		ofstream out(filename);
		for(i = 0; i < matrix.size(); ++i) {
			for(j = 0; j < matrix[0].size(); ++j)
				out<<matrix[i][j]<<" ";
			out<<endl;
		}
		out.close();
	}

} ///:~
