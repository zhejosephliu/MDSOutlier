/*
 * permute.cpp
 *
 *  Created on: Jul 24, 2010
 *      Author: Zhe Liu
 */

#include <iostream>
#include <vector>
#include <cmath>

#include "head.h"

using namespace std;

void Permute(Vec & array) {

	int i = 0;
	int length = array.size();
	int k = length;

	for(i = 0; i < length - 1; ++i) {
		double random = rand() / (double)(RAND_MAX);
		int index = (int)floor(random * k) + 1;
		int temp = array[index - 1];
		array[index - 1] = array[k - 1];
		array[k - 1] = temp;
		k--;
	}

} ///:~
