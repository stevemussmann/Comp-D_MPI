/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   DstatParent.h
 * Author: smuss
 *
 * Created on May 9, 2016, 8:03 AM
 */

#ifndef DSTATPARENT_H
#define DSTATPARENT_H

#include "Stats.h"
#include "fourtax.h"
#include "locusfile.h"

#include <fstream>
#include <iostream>

#include <math.h>

class DstatParent: public Stats{
    
    protected:
	void writeout(std::string *array, std::string output, int i, bool hetIgnore, bool hetInclude );
    public:
        DstatParent();
        DstatParent(const DstatParent& orig);
        virtual ~DstatParent();
        virtual void calcChiSqr() = 0;
	virtual void polyCalcChiSqr() = 0;
	virtual void calcDs(fourtax &dtest, unsigned int length, int ntaxa, locusfile &file, int my_rank) = 0;
	virtual void calcPolyDs(fourtax &dtest, int nloci) = 0;
	virtual void calcStats(unsigned int length) = 0;
	virtual void polyCalcStats(unsigned int length) = 0;
	virtual void polyCalcStats() = 0;
	virtual void calcStats() = 0;
	virtual void bootstrap(int mpiboot, int bootstrap, std::unordered_map <std::string,int> &indlist, 
			std::default_random_engine &generator, int ntaxa, locusfile &current, 
			std::vector<int> &keep, int my_rank, int i, int combs, std::string *indarray, 
			std::string output, bool hetIgnore, bool hetInclude) = 0;
	virtual void write(std::string *array, std::ofstream &outfile, int i, bool hetIgnore, bool hetInclude) = 0;
	//virtual void popCalcZ() = 0;

};

#endif /* DSTATPARENT_H */

