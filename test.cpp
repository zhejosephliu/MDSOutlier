/*
 * test.cpp
 *
 *  Created on: Jul 24, 2010
 *      Author: Zhe Liu
 */

#include <iostream>
#include <cmath>

#include "head.h"

using namespace std;

double PStest(double x11, double x12, double x21, double x22) {

	double xr1 = x11 + x12;
	double xr2 = x21 + x22;
	double xc1 = x11 + x21;
	double xc2 = x12 + x22;

  double xtotal = x11 + x12 + x21 + x22;

	double m11 = xr1 * xc1 / xtotal;
	double m12 = xr1 * xc2 / xtotal;
	double m21 = xr2 * xc1 / xtotal;
	double m22 = xr2 * xc2 / xtotal;

	if(m11 == 0 || m12 == 0 || m21 == 0 || m22 == 0) {
		return 0;
	}

	double ps_value = pow(x11 - m11, 2) / m11 +
                      pow(x12 - m12, 2) / m12 +
                      pow(x21 - m21, 2) / m21 +
                      pow(x22 - m22, 2) / m22;

	return ps_value;

} ///:~

double LRtest(double x11, double x12, double x21, double x22) {

	double xr1 = x11 + x12;
	double xr2 = x21 + x22;
	double xc1 = x11 + x21;
	double xc2 = x12 + x22;

	double xtotal = x11 + x12 + x21 + x22;

	double m11 = xr1 * xc1 / xtotal;
	double m12 = xr1 * xc2 / xtotal;
	double m21 = xr2 * xc1 / xtotal;
	double m22 = xr2 * xc2 / xtotal;

	if(m11 == 0 || m12 == 0 || m21 == 0 || m22 == 0) {
		return 0;
	}

	double c11;
	double c12;
	double c21;
	double c22;

	if(x11 == 0) {
		c11 = 0;
	}
	else {
		c11 = x11 * log(x11 / m11);
	}

	if(x12 == 0) {
		c12 = 0;
	}
	else {
		c12 = x12 * log(x12 / m12);
	}

	if(x21 == 0) {
		c21 = 0;
	}
	else {
		c21 = x21 * log(x21 / m21);
	}

	if(x22 == 0) {
		c22 = 0;
	}
	else {
		c22 = x22 * log(x22 / m22);
	}

	double lr_value = 2 * (c11 + c12 + c21 + c22);

	return lr_value;

} ///:~
