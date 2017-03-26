/*
 * mdsoutlier.cpp
 *
 *  Created on: Jul 24, 2010
 *      Author: Zhe Liu
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <ctime>
#include <cmath>

#include "head.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

using namespace std;

//: ------------------------------------------------------------------------------------------------------ //
//: mdsoutlier.cpp                                                                                         //
//: ------------------------------------------------------------------------------------------------------ //
int main(int argc, char** argv) {

	int numofcomp     =   atoi(argv[1]);           // 1st parameter: number of components in MDS
	int numofpermute  =   atoi(argv[2]);           // 2nd parameter: number of permutation times
	int numofremove   =   atoi(argv[3]);           // 3rd parameter: number of maximum outlyings

	cout<<endl;
	cout<<"number of components in MDS: "<<numofcomp;      cout<<endl;
	cout<<"number of permutation times: "<<numofpermute;   cout<<endl;
	cout<<"number of maximum outlyings: "<<numofremove;    cout<<endl;
	cout<<endl;

	// ----------------------------------------------------------------------------------------------------- //
	//                                    Read Data Set -- "Data.txt"                                        //
	// ----------------------------------------------------------------------------------------------------- //
	// Note: each row one snp, last row is the phenotype, coding 1 as case, and 0 as control                 //
	// Note: each column one individual, coding AA as 0, AB as 1, BB as 2, and missing as -1                 //
	// Note: there should not be a blank line at the end of the input file.                                  //
	// ----------------------------------------------------------------------------------------------------- //

	cout<<"Read Data Set."<<endl;

	ifstream in("Data.txt", ios::in);

	int snp = 0;                          // number of individuals
	int ids = 0;                          // number of SNPs

	Vec phenotype;                        // phenotype

	// ----------------------------------------------------------------------------------------------------- //

	for(string s; getline(in, s); ) {

		snp ++;

		istringstream sin(s);

		Vec line;
		for(double temp; sin>>temp; )
			line.push_back(temp);
		ids = line.size();

		int check;
		if(snp == 1) {
			check = ids;
			continue;
		}
		if(check != ids) {
			cout<<"ERROR!"<<endl;
			cout<<"The input dataset has false format."<<endl;
			exit(1);
		}

		phenotype = line;

	}

	// ----------------------------------------------------------------------------------------------------- //

	in.close();

	Vec tempphenotype = phenotype;       // temp phenotype

	snp --;                              // exclude phenotype

	cout<<endl;
	cout<<"number of individuals: "<<ids<<endl;
	cout<<"number of bio-markers: "<<snp<<endl;
	cout<<endl;

	// ----------------------------------------------------------------------------------------------------- //

	srand((unsigned)time(NULL));             // set seed for random number

	Mat ibs(ids, Vec(ids, 0));               // IBS matrix: ids * ids

	Mat coordinate(ids, Vec(numofcomp, 0));    // coordinate matrix: individuals * numofcomp


	// record FINAL P-VALUES in observed data and permutation samples
	Vec observe_pvalue(numofremove + 1,0);
	Mat permute_pvalue(numofpermute, Vec((numofremove + 1), 0));


	// record CHI-SQUARES in "original data" and "best-removal data"
	Vec original_chisq(snp, 0);
	Vec best_chisq(snp, 0);


	// record PVALUES in "original data" and "best-removal data"
	Vec original_pvalue(snp, 0);
	Vec best_pvalue(snp, 0);


	int num_best;        // the order of the best remove
	double value_best;   // the minimal observed p-value


	// record the ID of individuals
	Vec ids_remove;
	for(int r = 0; r < ids; ++r) {
		ids_remove.push_back(r + 1);
	}


	// record the genotypes of removals
	Mat removegenotype(snp, Vec(numofremove, 0));


	// record the order of removed individuals
	Vec removeindividual;


	// ----------------------------------------------------------------------------------------------------- //
	//                                               Calculation                                             //
	// ----------------------------------------------------------------------------------------------------- //

	cout<<"Calculation."<<endl;

	cout<<endl;
	cout<<"Compute the IBS Dissimilarity Matrix."<<endl;


	// IBS analysis
	IBSmatrix(snp, ids, ibs);
	Printmatrix(ibs, "IBSmatrix.txt");


	// The following is to import the IBS matrix that has already been computed.
	// Be sure to make the IBS calculation above into comments before.


	/*---------------------------------------------
	ifstream inputIBS("IBSmatrix.txt", ios::in);
	int cnum = 0;
	for(string s; getline(inputIBS, s); ) {
      Vec line;
      istringstream sin(s);
      for(double temp; sin>>temp; )
      line.push_back(temp);
      if(ibs[cnum].size() != line.size()) {
        cout<<"Error: Not Equal Size"<<endl;
		exit(1);
      }
      ibs[cnum] = line;
      cnum ++;
	}
	inputIBS.close();
	---------------------------------------------*/

	cout<<endl;
	cout<<"Observed Data Analysis: Allelic PS Test."<<endl;


	// ----------------------------------------------------------------------------------------------------- //


	Mat table(snp, Vec(4, 0));      // 2-by-2 table for each SNP
	Vec chisq(snp, 0);              // chi-square value for each SNP

	ifstream input("Data.txt", ios::in);

	int indicator = -1;
	for(string s; getline(input, s); ) {
		indicator ++;
		if(indicator == snp) {
			break;
		}
		Vec temptable(6, 0);
		istringstream sin(s);
		int num = -1;
		for(double temp; sin>>temp; ) {
			num ++;
			if(phenotype[num] == 1) {
				if(temp != -1)
					temptable[temp] ++;
			}
			else {
				if(phenotype[num] != 0) {
					cout<<"ERROR!"<<endl;
					cout<<"Phenotypes are not coded with 1 and 0"<<endl;
					exit(1);
				}
				if(temp != -1)
					temptable[temp+3] ++;
			}
		}
		double x11 = 2 * temptable[0] + temptable[1];
		double x12 = 2 * temptable[2] + temptable[1];
		double x21 = 2 * temptable[3] + temptable[4];
		double x22 = 2 * temptable[5] + temptable[4];
		table[indicator][0] = x11;
		table[indicator][1] = x12;
		table[indicator][2] = x21;
		table[indicator][3] = x22;
		chisq[indicator] = PStest(x11, x12, x21, x22);
	}

	input.close();

	// ----------------------------------------------------------------------------------------------------- //


	// print original chi-squares and p-values
	original_chisq = chisq;
	for(int c = 0; c < snp; ++c) {
		original_pvalue[c] = 1 - gsl_cdf_chisq_P(original_chisq[c], 1);
	}

	// Printvector(original_chisq, "original_chisq.txt");
	Printvector(original_pvalue, "original_pvalue.txt");

	// minimum p-value
	observe_pvalue[0] = *min_element(original_pvalue.begin(), original_pvalue.end());

	cout<<endl;
	cout<<"Observed Data Minimal P-value: "<<observe_pvalue[0]<<"."<<endl;


	// ----------------------------------------------------------------------------------------------------- //
	//                                                 Remove                                                //
	// ----------------------------------------------------------------------------------------------------- //

	cout<<endl;
	cout<<"Remove Individuals"<<endl<<endl;

	best_chisq   =  original_chisq;

	num_best     =  0;

	value_best   =  observe_pvalue[0];

	// ----------------------------------------------------------------------------------------------------- //

	for(int n=1; n <= numofremove; ++n) {

		ifstream input("Data.txt", ios::in);

		// MDS analysis
		MDS(ibs, coordinate, numofcomp);
		if(n == 1) {
			Printmatrix(coordinate, "MDSmatrix.txt");
		}

		// which one to remove
		int remove;
		double maxdistance;
		for(int z = 0; z < ids - n + 1; ++z) {
			double sumdistance = 0;
			for(int x = 0; x < numofcomp; ++x)
				sumdistance += pow(coordinate[z][x], 2);
			sumdistance = sqrt(sumdistance);
			if(z == 0) {
				remove = 1;
				maxdistance = sumdistance;
			}
			else {
				if(sumdistance > maxdistance) {
					remove = z + 1;
					maxdistance = sumdistance;
				}
			}
		}

		// remove: from "1" to "ids"
		removeindividual.push_back(ids_remove[remove - 1]);

		// output
		cout<<n<<": remove individual "<<ids_remove[remove - 1]<<"."<<endl;

		// --------------------------------------------------------------------------------------------------- //

		// Data Analysis
		int indicator = -1;
		for(string s; getline(input, s); ) {
			indicator ++;
			if(indicator == snp) {
				break;
			}
			istringstream sin(s);
			int num = -1;
			for(double temp; sin>>temp; ) {
				num ++;
				if(num != ids_remove[remove - 1] - 1) {
					continue;
				}
				removegenotype[indicator][n - 1] = temp;   // this is used for permutation test
				if(phenotype[num] == 1) {
					if(temp == 0)
						table[indicator][0] -= 2;
					if(temp == 1) {
						table[indicator][0] -= 1;
						table[indicator][1] -= 1;
					}
					if(temp == 2)
						table[indicator][1] -= 2;
				}
				else {
					if(phenotype[num] != 0) {
						cout<<"ERROR!"<<endl;
						cout<<"Phenotypes are not coded with 1 and 0"<<endl;
						exit(1);
					}
					if(temp == 0)
						table[indicator][2] -= 2;
					if(temp == 1) {
						table[indicator][2] -= 1;
						table[indicator][3] -= 1;
					}
					if(temp == 2)
						table[indicator][3] -= 2;
				}
				break;
			}
			chisq[indicator] = PStest(table[indicator][0], table[indicator][1], table[indicator][2], table[indicator][3]);
		}

		input.close();

		// --------------------------------------------------------------------------------------------------- //

		// maximal remove-chi-square
		double rmaxchisq = *max_element(chisq.begin(), chisq.end());

		// minimal remove-p-value
		observe_pvalue[n] = 1 - gsl_cdf_chisq_P(rmaxchisq, 1);

		// output
		cout<<n<<": p-value "<<observe_pvalue[n]<<"."<<endl<<endl;

		// select the best removal and its p-values
		if(observe_pvalue[n] < value_best) {
			num_best = n;
			value_best = observe_pvalue[n];
			best_chisq = chisq;
		}

		// --------------------------------------------------------------------------------------------------- //

		// update the remove array
		ids_remove.erase(ids_remove.begin() + remove - 1);

		// update the IBS matrix: remove the (remove)th row and column
		ibs.erase(ibs.begin() + remove - 1);
		for(unsigned int ri = 0; ri < ibs.size(); ++ri)
			ibs[ri].erase(ibs[ri].begin() + remove - 1);

		// update coordinate matrix: remove the (remove)th row
		coordinate.erase(coordinate.begin() + remove - 1);

	}

	for(int c = 0; c < snp; ++c) {
		best_pvalue[c] = 1 - gsl_cdf_chisq_P(best_chisq[c], 1);
	}

	Printvector(observe_pvalue, "ObservePvalueTrend.txt");

	Printvector(best_chisq, "best_chisq.txt");

	Printvector(best_pvalue, "best_pvalue.txt");

	Printvector(removeindividual, "RemovedIDs.txt");

	Printmatrix(removegenotype, "remvoegenotype.txt");

	Printmatrix(coordinate, "BestMDSmatrix.txt");

	// ----------------------------------------------------------------------------------------------------- //
	//                                                 Permutation                                           //
	// ----------------------------------------------------------------------------------------------------- //


	for(int p = 0; p < numofpermute; ++p) {

		cout<<"Permutation "<<(p + 1)<<" with minimal p-value: ";

		Mat ptable(snp, Vec(4, 0));
		Vec pchisq(snp, 0);

		Vec pphenotype = tempphenotype;
		Permute(pphenotype);

		ifstream input("Data.txt", ios::in);

		int indicator = -1;
		for(string s; getline(input, s); ) {
			indicator ++;
			if(indicator == snp) {
				break;
			}
			Vec ptemptable(6, 0);
			istringstream sin(s);
			int num = -1;
			for(double temp; sin>>temp; ) {
				num ++;
				if(pphenotype[num] == 1) {
					if(temp != -1)
						ptemptable[temp] ++;
				}
				else {
					if(pphenotype[num] != 0) {
						cout<<"ERROR!"<<endl;
						cout<<"Phenotypes are not coded with 1 and 0"<<endl;
						exit(1);
					}
					if(temp != -1)
						ptemptable[temp + 3] ++;
				}
			}
			double x11 = 2 * ptemptable[0] + ptemptable[1];
			double x12 = 2 * ptemptable[2] + ptemptable[1];
			double x21 = 2 * ptemptable[3] + ptemptable[4];
			double x22 = 2 * ptemptable[5] + ptemptable[4];
			ptable[indicator][0] = x11;
			ptable[indicator][1] = x12;
			ptable[indicator][2] = x21;
			ptable[indicator][3] = x22;
			pchisq[indicator] = PStest(x11, x12, x21, x22);
		}

		input.close();

		double pmaxchisq = *max_element(pchisq.begin(), pchisq.end());

		permute_pvalue[p][0] = 1 - gsl_cdf_chisq_P(pmaxchisq, 1);

		for(int pn = 1; pn <= numofremove; ++pn) {
			if(pphenotype[(removeindividual[pn - 1]) - 1] == 1) {
				for(int ps = 0; ps < snp; ++ps) {
					if(removegenotype[ps][pn - 1] == 0)
						ptable[ps][0] -= 2;
					if(removegenotype[ps][pn - 1] == 1) {
						ptable[ps][0] -= 1;
						ptable[ps][1] -= 1;
					}
					if(removegenotype[ps][pn - 1] == 2)
						ptable[ps][1] -= 2;
					pchisq[ps] = PStest(ptable[ps][0], ptable[ps][1], ptable[ps][2], ptable[ps][3]);
				}
			}
			else {
				if(pphenotype[(removeindividual[pn - 1]) - 1] != 0) {
					cout<<"ERROR!"<<endl;
					cout<<"Phenotypes are not coded with 1 and 0"<<endl;
					exit(1);
				}
				for(int ps = 0; ps < snp; ++ps) {
					if(removegenotype[ps][pn - 1] == 0)
						ptable[ps][2] -= 2;
					if(removegenotype[ps][pn - 1] == 1) {
						ptable[ps][2] -= 1;
						ptable[ps][3] -= 1;
					}
					if(removegenotype[ps][pn - 1] == 2)
						ptable[ps][3] -= 2;
					pchisq[ps] = PStest(ptable[ps][0], ptable[ps][1], ptable[ps][2], ptable[ps][3]);
				}
			}
			pmaxchisq = *max_element(pchisq.begin(), pchisq.end());
			permute_pvalue[p][pn] = 1 - gsl_cdf_chisq_P(pmaxchisq, 1);
		}

		cout<<(*min_element(permute_pvalue[p].begin(), permute_pvalue[p].end()))<<endl;

	}

	Printmatrix(permute_pvalue, "permute_pvalue.txt");

	// ----------------------------------------------------------------------------------------------------- //
	//                                                  Finish                                               //
	// ----------------------------------------------------------------------------------------------------- //

	// calculate the overall pvalue
	double pvalue = 0;
	double originpvalue = 0;

	Vec permute_best(numofpermute, 0);

	// calculate the empirical pvalue
	Vec original_empirical(snp, 0);
	Vec best_empirical(snp, 0);

	for(int k = 0; k < numofpermute; ++k) {

		// MIN_PVALUE
		// permute_best[k] = *min_element(permute_pvalue[k].begin(), permute_pvalue[k].end());

		// MIN_INTER_PVALUE
		unsigned int f;
		for(f = 1; f < permute_pvalue[k].size(); ++f) {
			if(permute_pvalue[k][f] > permute_pvalue[k][f - 1]) {
				continue;
			}
			else {
				break;
			}
		}
		if(f == permute_pvalue[k].size()) {
			permute_best[k] = permute_pvalue[k][f - 1];
		}
		else {
			permute_best[k] = permute_pvalue[k][f];
			for(unsigned int w = f; w < permute_pvalue[k].size(); ++w) {
				if(permute_pvalue[k][w] < permute_best[k]) {
					permute_best[k] = permute_pvalue[k][w];
				}
			}
		}

		for(int y = 0; y < snp; ++y) {
			if(original_pvalue[y] >= permute_pvalue[k][0])
				original_empirical[y] ++;
			if(best_pvalue[y] >= permute_best[k])
				best_empirical[y] ++;
		}

		if(value_best >= permute_best[k])
			pvalue ++;

		if(observe_pvalue[0] >= permute_pvalue[k][0])
			originpvalue ++;

	}

	for(int z = 0; z < snp; ++z) {
		original_empirical[z] /= numofpermute;
		best_empirical[z] /= numofpermute;
	}

	Printvector(original_empirical, "original_empirical.txt");
	Printvector(best_empirical, "best_empirical.txt");

	originpvalue /= numofpermute;
	pvalue /= numofpermute;

	// output result
	cout<<endl;

	cout<<"The order of best removal is: "<<num_best<<endl;

	cout<<endl;

	cout<<"The original p-value is: "<<originpvalue<<endl;

	cout<<endl;

	cout<<"The overall p-value is: "<<pvalue<<endl;

	return 0;

} ///:~
