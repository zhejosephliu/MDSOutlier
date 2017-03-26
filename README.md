# mdsOutlier

MDSOutlier is a whole-genome case-control association analysis program, designed for systematic removal of outliers to reduce heterogeneity.

For more information please contact zhe.liu.uchicago@gmail.com

#### Description

In human case-control association studies, population heterogeneity is often present and can lead to increased false-positive results. Various methods have been proposed and are in current use to remedy this situation. We assume that heterogeneity is due to a relatively small number of individuals whose allele frequencies differ from those of the remainder of the sample. For this situation, we propose a new method of handling heterogeneity by removing outliers in a controlled manner. In a coordinate system of the c largest principal components in multidimensional scaling (MDS), we systematically remove one after another of the most extreme outlying individuals and each time recompute the largest association test statistic. The smallest p-value obtained within M removals serves as our test statistic whose significance level is assessed in randomization samples. In power simulations of our method and three methods in current use, averaged over several different scenarios, the best method turned out to be logistic regression analysis (based on all individuals) with MDS components as covariates. Our proposed method ranked closely behind logistic regression analysis with MDS components but ahead of other commonly used approaches. In analyses of real datasets our method performed best.

#### Installation 

Only for Linux/Unix and Mac OS X platforms.

GNU Scientific Library (GSL) should be installed before the compilation.

    $g++ -Wall -I/usr/local/include -c Main.cpp IBSmatrix.cpp MDS.cpp Print.cpp Permute.cpp Test.cpp
    $g++ -static Main.o IBSmatrix.o MDS.o Print.o Permute.o Test.o -lgsl -lgslcblas -lm -o MDSOutlier

#### Usage

    $./MDSOutlier component permute removal

    -- MDSOutlier:  the executable file
    -- component:   number of principal components, e.g. 4
    -- permute:     number of permutations, e.g. 5000
    -- removal:     number of maximal outliers, e.g. 20
    
    -- Data.txt:    re-format the data according to the following instruction,
                    re-name it to "Data.txt", and put it into the same directory with the program);
                    each row (except the last row) represents genotype information of one snp, 
                    while the last row represents phenotype information (case = 1, control = 0);
                    each column reperesents one individual, coding AA as 0, AB as 1, BB as 2, and
                    missing as -1, where the allele B is the minor allele; e.g. 6 individuals and 4 SNPs:
                    2 1 0 1 0 1
                    1 2 2 0 1 1
                    2 1 1 1 0 2
                    0 2 0 2 1 2
                    1 1 1 0 0 0

#### References

Y Shen, Z Liu, and J Ott. Systematic removal of outliers to reduce heterogeneity in case-control association studies. Human Heredity, 70(4):227-231, 2010.
