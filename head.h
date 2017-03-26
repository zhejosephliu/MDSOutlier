/*
 * head.h
 *
 *  Created on: Jul 24, 2010
 *      Author: Zhe Liu
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <list>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <cstring>

using namespace std;

typedef vector<double> Vec;
typedef vector<vector<double> > Mat;

//: print.cpp
void Printvector(const Vec &, const char *); ///:~
void Printmatrix(const Mat &, const char *); ///:~

//: ibsmatrix.cpp
void IBSmatrix(const int, const int, Mat &); ///:~

//: test.cpp
double PStest(double, double, double, double); ///:~
double LRtest(double, double, double, double); ///:~

//: permute.cpp
void Permute(Vec &); ///:~

//: mds.cpp
void MDS(Mat &, Mat &, int); ///:~
