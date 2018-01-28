/* 
 * File:   Dstat.h
 * Author: Steve
 *
 * Created on July 26, 2015, 11:48 AM
 */

#ifndef DSTAT_H
#define	DSTAT_H

#include "fourtax.h"
#include "DstatParent.h"
#include "Stats.h"
#include "locusfile.h"

#include <vector>
#include <unordered_map>
#include <random>
#include <fstream>
#include <iostream>

class Dstat: public DstatParent{
public:
	Dstat();
	//integer based method of calculating D-statistics
	void calcDstat(fourtax &dtest, unsigned int length, int ntaxa);
	
	//double-based method of calculating D-statistics
	void polyCalcDstat(fourtax &dtest, unsigned int length, int ntaxa);
	
	double getD();
	int getABBA();
	int getBABA();
	double getChiSqr();
	Dstat(const Dstat& orig);
	virtual ~Dstat();
	
	//overridden function sfrom parent class
	void calcChiSqr() override;
	void polyCalcChiSqr() override;
	void calcDs(fourtax &dtest, unsigned int length, int ntaxa, locusfile &file, int my_rank) override;
	void calcPolyDs(fourtax &dtest, int nloci);
	void calcStats(unsigned int length) override;
	void polyCalcStats(unsigned int length) override;
	void calcStats() override;
	void polyCalcStats() override;
	void bootstrap(int mpiboot, int bootstrap, std::unordered_map <std::string,int> &indlist, 
		   std::default_random_engine &generator, int ntaxa, locusfile &current, 
		   std::vector<int> &keep, int my_rank, int i, int combs, std::string *indarray, 
		   std::string output, bool hetIgnore, bool hetInclude) override;
	void write(std::string *array, std::ofstream &outfile, int i, bool hetIgnore, bool hetInclude) override;
	//void popCalcZ() override;
private:
	double D;
	int ABBA;
	int BABA;
	unsigned int numloci;
	double polyABBA;
	double polyBABA;
	double chisq;
	double avg;
	double sd;
	double Z;
	double Zpval;
	double Xpval;
	void bootproc(int bootstrap, std::vector<int> &keep, double *bootD, locusfile &file, 
		  std::unordered_map <std::string,int> &indlist, std::default_random_engine &generator, 
		  int ntaxa, bool hetIgnore, bool hetInclude, int my_rank);
};

#endif	/* DSTAT_H */

