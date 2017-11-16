/* 
 * File:   partD.h
 * Author: Steve
 *
 * Created on July 23, 2015, 8:22 PM
 */

#ifndef PARTD_H
#define	PARTD_H

#include "fourtax.h"
#include "Stats.h"
#include "locusfile.h"

#include <unordered_map>
#include <string>
#include "DstatParent.h"
#include <vector>
#include <random>
#include <fstream>
#include <iostream>

class partD: public DstatParent {
public:
    partD();
    //for integer-based calculations of D-statistics
    void calcD12();
    void calcD1();
    void calcD2();
    
    //for double-based calculations of D-statistics
    void polyCalcD12();
    void polyCalcD1();
    void polyCalcD2();
    
    double getD1();
    double getD2();
    double getD12();
    double getChiSqrD1();
    double getChiSqrD2();
    double getChiSqrD12();
    void calcPartD(fourtax &dtest, unsigned int length, int ntaxa);
    int getPatternInt(std::string pattern);
    partD(const partD& orig);
    virtual ~partD();
    
    //overridden functions from parent class
    void calcChiSqr() override;
    void polyCalcChiSqr() override;
    void calcDs(fourtax &dtest, unsigned int length, int ntaxa, locusfile &file, int my_rank) override;
    void calcPolyDs(fourtax &dtest, int nloci);
    void calcStats(unsigned int length) override;
    void polyCalcStats(unsigned int length) override;
    void calcStats() override;
    void polyCalcStats() override;
    void bootstrap(int mpiboot, int bootstrap, std::unordered_map <std::string,int> &indlist, 
		std::default_random_engine &generator, int ntaxa, locusfile &current, std::vector<int> &keep, 
		int my_rank, int i, int combs, std::string *indarray, std::string output, bool hetIgnore, bool hetInclude) override;
    void write(std::string *array, std::ofstream &outfile, int i, bool hetIgnore, bool hetInclude) override;
    //void popCalcZ() override;
    
private:
    double D12;
    double D1;
    double D2;
    double xsqrD1;
    double xsqrD2;
    double xsqrD12;
    double avgD1;
    double avgD2;
    double avgD12;
    double sdD1;
    double sdD2;
    double sdD12;
    double ZD1;
    double ZD2;
    double ZD12;
    double ZD1pval;
    double ZD2pval;
    double ZD12pval;
    double XpvalD1;
    double XpvalD2;
    double XpvalD12;
    std::unordered_map <std::string,int> patterns;
    std::unordered_map<std::string,double> polyPatterns;
    void bootproc(int bootstrap, std::vector<int> &keep, double *bootD1, double *bootD2, double *bootD12, 
		  locusfile &file, std::unordered_map <std::string,int> &indlist, 
		  std::default_random_engine &generator, int ntaxa, bool hetIgnore, bool hetInclude, int my_rank);
};

#endif	/* PARTD_H */

