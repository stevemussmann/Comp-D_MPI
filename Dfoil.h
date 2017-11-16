/* 
 * File:   Dfoil.h
 * Author: Steve
 *
 * Created on July 26, 2015, 7:24 AM
 */

#ifndef DFOIL_H
#define	DFOIL_H

#include "fourtax.h"
#include "Stats.h"
#include "locusfile.h"

#include <unordered_map>
#include <string>

#include "DstatParent.h"

#include <boost/math/distributions/chi_squared.hpp>

class Dfoil: public DstatParent {
public:
    Dfoil();
    //integer based methods for calculating D-statistics
    void calcDFO();
    void calcDIL();
    void calcDFI();
    void calcDOL();
    
    //double-based methods for calculating D-statistics
    void polyCalcDFO();
    void polyCalcDIL();
    void polyCalcDFI();
    void polyCalcDOL();
    
    double getDFO();
    double getDIL();
    double getDFI();
    double getDOL();
    double getChiSqrDFO();
    double getChiSqrDIL();
    double getChiSqrDFI();
    double getChiSqrDOL();
    void calcDfoil(fourtax &dtest, unsigned int length, int ntaxa);
    int getPatternInt(std::string pattern);
    Dfoil(const Dfoil& orig);
    virtual ~Dfoil();
    
    //overridden functions from parent class
    void calcChiSqr() override;
    void polyCalcChiSqr() override;
    void calcDs(fourtax &dtest, unsigned int length, int ntaxa, locusfile &file, int my_rank) override;
    void calcPolyDs(fourtax &dtest, int nloci);
    void calcStats(unsigned int length) override;
    void polyCalcStats(unsigned int length) override;
    void polyCalcStats() override;
    void calcStats() override;
    void bootstrap(int mpiboot, int bootstrap, std::unordered_map <std::string,int> &indlist, 
		   std::default_random_engine &generator, int ntaxa, locusfile &current, std::vector<int> &keep, 
		   int my_rank, int i, int combs, std::string *indarray, std::string output, 
		   bool hetIgnore, bool hetInclude) override;
    void write(std::string *array, std::ofstream &outfile, int i, bool hetIgnore, bool hetInclude) override;
    //void popCalcZ() override;
    
private:
    double DFO;
    double DIL;
    double DFI;
    double DOL;
    double xsqrDFO;
    double xsqrDIL;
    double xsqrDFI;
    double xsqrDOL;
    double avgDFO;
    double avgDIL;
    double avgDFI;
    double avgDOL;
    double sdDFO;
    double sdDIL;
    double sdDFI;
    double sdDOL;
    double ZDFO;
    double ZDIL;
    double ZDFI;
    double ZDOL;
    double ZDFOpval;
    double ZDILpval;
    double ZDFIpval;
    double ZDOLpval;
    double XpvalDFO;
    double XpvalDIL;
    double XpvalDFI;
    double XpvalDOL;
    std::unordered_map <std::string,int> patterns;
    std::unordered_map<std::string,double> polyPatterns;
    void bootproc(int bootstrap, std::vector<int> &keep, double *bootDFO, double *bootDIL, double *bootDFI, 
		  double *bootDOL, locusfile &file, std::unordered_map <std::string,int> &indlist, 
		  std::default_random_engine &generator, int ntaxa, bool hetIgnore, bool hetInclude, int my_rank);
};

#endif	/* DFOIL_H */

