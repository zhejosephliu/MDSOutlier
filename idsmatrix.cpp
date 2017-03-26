/*
 * ibsmatrix.cpp
 *
 *  Created on: Jul 24, 2010
 *      Author: Zhe Liu
 */

#include <iostream>
#include <vector>
#include <cmath>

#include "head.h"

using namespace std;

void IBSmatrix(const int snp, const int ids, Mat & ibs) {

	unsigned int i, j;

	int k = 0;

	int length = ibs.size();

	Mat nonmissing = ibs;

	if(length != ids) {
		cout<<"ERROR!"<<endl;
		cout<<"Numbers of row and column are not equal."<<endl;
		exit(1);
	}

	ifstream in("Data.txt", ios::in);

	for(string s; getline(in, s); ) {
		k++;
		if(k == snp + 1) {
			break;
		}
		Vec line;
		istringstream sin(s);
		for(double temp; sin>>temp; )
			line.push_back(temp);
		for(i = 0; i < ibs.size(); ++i) {
			for(j = i + 1; j < ibs.size(); ++j) {
				if(line[i] == -1 || line[j] == -1) {
					continue;
				}
				nonmissing[i][j] ++;
				double compare = abs(line[i] - line[j]);
				if(compare == 2) ibs[i][j] += 0;
				if(compare == 1) ibs[i][j] += 0.5;
				if(compare == 0) ibs[i][j] += 1;
			}
		}
	}

	in.close();

	for(i = 0; i < ibs.size(); ++i) {
		for(j = i + 1; j < ibs.size(); ++j) {
			ibs[i][j] = 1 - ibs[i][j] / nonmissing[i][j];
			ibs[j][i] = ibs[i][j];
		}
		ibs[i][i] = 0;
	}

} ///:~
