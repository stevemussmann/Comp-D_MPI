/* 
 * File:   Dfoil.cpp
 * Author: Steve
 * 
 * Created on July 26, 2015, 7:24 AM
 */

#include "Dfoil.h"
#include "fourtax.h"
#include "Stats.h"
#include "locusfile.h"

#include <unordered_map>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "mpi.h"

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>

Dfoil::Dfoil() {
	DFO=0.0;
	DIL=0.0;
	DFI=0.0;
	DOL=0.0;
	patterns["AAABA"] = 0;
	patterns["AABAA"] = 0;
	patterns["ABAAA"] = 0;
	patterns["ABABA"] = 0;
	patterns["ABBAA"] = 0;
	patterns["ABBBA"] = 0;
	patterns["BAAAA"] = 0;
	patterns["BAABA"] = 0;
	patterns["BABAA"] = 0;
	patterns["BABBA"] = 0;
	patterns["BBABA"] = 0;
	patterns["BBBAA"] = 0;
	
	//for calculating dstat incorporating heterozygous individuals
	polyPatterns["AAABA"] = 0.0;
	polyPatterns["AABAA"] = 0.0;
	polyPatterns["ABAAA"] = 0.0;
	polyPatterns["ABABA"] = 0.0;
	polyPatterns["ABBAA"] = 0.0;
	polyPatterns["ABBBA"] = 0.0;
	polyPatterns["BAAAA"] = 0.0;
	polyPatterns["BAABA"] = 0.0;
	polyPatterns["BABAA"] = 0.0;
	polyPatterns["BABBA"] = 0.0;
	polyPatterns["BBABA"] = 0.0;
	polyPatterns["BBBAA"] = 0.0;
	
	avgDFO = 0.0;
	avgDIL = 0.0;
	avgDFI = 0.0;
	avgDOL = 0.0;

	sdDFO = 0.0;
	sdDIL = 0.0;
	sdDFI = 0.0;
	sdDOL = 0.0;

	ZDFO = 0.0;
	ZDIL = 0.0;
	ZDFI = 0.0;
	ZDOL = 0.0;

	ZDFOpval = 0.0;
	ZDILpval = 0.0;
	ZDFIpval = 0.0;
	ZDOLpval = 0.0;

	XpvalDFO = 0.0;
	XpvalDIL = 0.0;
	XpvalDFI = 0.0;
	XpvalDOL = 0.0;

	numloci=0;	

}

Dfoil::Dfoil(const Dfoil& orig) {
}

Dfoil::~Dfoil() {
}

void Dfoil::calcDFO()
{
	int denom = (patterns["BABAA"]+patterns["BBBAA"]+patterns["ABABA"]+patterns["AAABA"])+(patterns["BAABA"]+patterns["BBABA"]+patterns["ABBAA"]+patterns["AABAA"]);
	if( denom == 0)
	{
		DFO=0.0;
	}
	else
	{
		int num=(patterns["BABAA"]+patterns["BBBAA"]+patterns["ABABA"]+patterns["AAABA"])-(patterns["BAABA"]+patterns["BBABA"]+patterns["ABBAA"]+patterns["AABAA"]);
		DFO=(double)num/(double)denom;
	}
}

void Dfoil::calcDIL()
{
	int denom = (patterns["ABBAA"]+patterns["BBBAA"]+patterns["BAABA"]+patterns["AAABA"])+(patterns["ABABA"]+patterns["BBABA"]+patterns["BABAA"]+patterns["AABAA"]);
	if(denom == 0)
	{
		DIL=0.0;
	}
	else
	{
		int num = (patterns["ABBAA"]+patterns["BBBAA"]+patterns["BAABA"]+patterns["AAABA"])-(patterns["ABABA"]+patterns["BBABA"]+patterns["BABAA"]+patterns["AABAA"]);
		DIL=(double)num/(double)denom;
	}
}

void Dfoil::calcDFI()
{
	int denom = (patterns["BABAA"]+patterns["BABBA"]+patterns["ABABA"]+patterns["ABAAA"])+(patterns["ABBAA"]+patterns["ABBBA"]+patterns["BAABA"]+patterns["BAAAA"]);
	if(denom == 0)
	{
		DFI=0.0;
	}
	else
	{
		int num=(patterns["BABAA"]+patterns["BABBA"]+patterns["ABABA"]+patterns["ABAAA"])-(patterns["ABBAA"]+patterns["ABBBA"]+patterns["BAABA"]+patterns["BAAAA"]);
		DFI=(double)num/(double)denom;
	}
}

void Dfoil::calcDOL()
{
	int denom = (patterns["BAABA"]+patterns["BABBA"]+patterns["ABBAA"]+patterns["ABAAA"])+(patterns["ABABA"]+patterns["ABBBA"]+patterns["BABAA"]+patterns["BAAAA"]);
	if(denom == 0)
	{
		DOL=0.0;
	}
	else
	{
		int num=(patterns["BAABA"]+patterns["BABBA"]+patterns["ABBAA"]+patterns["ABAAA"])-(patterns["ABABA"]+patterns["ABBBA"]+patterns["BABAA"]+patterns["BAAAA"]);
		DOL=(double)num/(double)denom;
	}
}

void Dfoil::polyCalcDFO()
{
	double denom = (polyPatterns["BABAA"]+polyPatterns["BBBAA"]+polyPatterns["ABABA"]+polyPatterns["AAABA"])+(polyPatterns["BAABA"]+polyPatterns["BBABA"]+polyPatterns["ABBAA"]+polyPatterns["AABAA"]);
	if( (int)denom == 0)
	{
		DFO=0.0;
	}
	else
	{
		double num=(polyPatterns["BABAA"]+polyPatterns["BBBAA"]+polyPatterns["ABABA"]+polyPatterns["AAABA"])-(polyPatterns["BAABA"]+polyPatterns["BBABA"]+polyPatterns["ABBAA"]+polyPatterns["AABAA"]);
		DFO=num/denom;
	}
}

void Dfoil::polyCalcDIL()
{
	double denom = (polyPatterns["ABBAA"]+polyPatterns["BBBAA"]+polyPatterns["BAABA"]+polyPatterns["AAABA"])+(polyPatterns["ABABA"]+polyPatterns["BBABA"]+polyPatterns["BABAA"]+polyPatterns["AABAA"]);
	if( (int)denom == 0)
	{
		DIL=0.0;
	}
	else
	{
		double num = (polyPatterns["ABBAA"]+polyPatterns["BBBAA"]+polyPatterns["BAABA"]+polyPatterns["AAABA"])-(polyPatterns["ABABA"]+polyPatterns["BBABA"]+polyPatterns["BABAA"]+polyPatterns["AABAA"]);
		DIL=num/denom;
	}
}

void Dfoil::polyCalcDFI()
{
	double denom = (polyPatterns["BABAA"]+polyPatterns["BABBA"]+polyPatterns["ABABA"]+polyPatterns["ABAAA"])+(polyPatterns["ABBAA"]+polyPatterns["ABBBA"]+polyPatterns["BAABA"]+polyPatterns["BAAAA"]);
	if( (int)denom == 0)
	{
		DFI=0.0;
	}
	else
	{
		double num=(polyPatterns["BABAA"]+polyPatterns["BABBA"]+polyPatterns["ABABA"]+polyPatterns["ABAAA"])-(polyPatterns["ABBAA"]+polyPatterns["ABBBA"]+polyPatterns["BAABA"]+polyPatterns["BAAAA"]);
		DFI=num/denom;
	}
}

void Dfoil::polyCalcDOL()
{
	double denom = (polyPatterns["BAABA"]+polyPatterns["BABBA"]+polyPatterns["ABBAA"]+polyPatterns["ABAAA"])+(polyPatterns["ABABA"]+polyPatterns["ABBBA"]+polyPatterns["BABAA"]+polyPatterns["BAAAA"]);
	if( (int)denom == 0)
	{
		DOL=0.0;
	}
	else
	{
		double num=(polyPatterns["BAABA"]+polyPatterns["BABBA"]+polyPatterns["ABBAA"]+polyPatterns["ABAAA"])-(polyPatterns["ABABA"]+polyPatterns["ABBBA"]+polyPatterns["BABAA"]+polyPatterns["BAAAA"]);
		DOL=num/denom;
	}
}

void Dfoil::calcDs(fourtax &dtest, unsigned int length, int ntaxa, locusfile &file, int my_rank)
{
	for(unsigned int i=0; i<length; i++)
	{
		dtest.calculatePattern(i, ntaxa, file, my_rank);
		patterns[dtest.getPattern(i)]++;
	}
}

void Dfoil::calcPolyDs(fourtax &dtest, int nloci)
{
	for(int i=0; i<nloci; i++)
	{
  
		double P1=dtest.getFreq(i, 4);
		double P2=dtest.getFreq(i, 3);
		double P3a=dtest.getFreq(i, 2);
		double P3b=dtest.getFreq(i, 1);
		double O=dtest.getFreq(i, 0);
	
	
		polyPatterns["ABBBA"]+=((1.0-P1)*P2*P3a*P3b*(1.0-O));
		polyPatterns["BABBA"]+=(P1*(1.0-P2)*P3a*P3b*(1.0-O));
		polyPatterns["ABBAA"]+=((1.0-P1)*P2*P3a*(1.0-P3b)*(1.0-O));
		polyPatterns["BABAA"]+=(P1*(1.0-P2)*P3a*(1.0-P3b)*(1.0-O));
		polyPatterns["ABABA"]+=((1.0-P1)*P2*(1.0-P3a)*P3b*(1.0-O));
		polyPatterns["BAABA"]+=(P1*(1.0-P2)*(1.0-P3a)*P3b*(1.0-O));
	
		polyPatterns["BBBAA"]+=(P1*P2*P3a*(1.0-P3b)*(1.0-O));
		polyPatterns["BBABA"]+=(P1*P2*(1.0-P3a)*P3b*(1.0-O));
	
		polyPatterns["AAABA"]+=((1.0-P1)*(1.0-P2)*(1.0-P3a)*P3b*(1.0-O));
		polyPatterns["AABAA"]+=((1.0-P1)*(1.0-P2)*P3a*(1.0-P3b)*(1.0-O));
		polyPatterns["ABAAA"]+=((1.0-P1)*P2*(1.0-P3a)*(1.0-P3b)*(1.0-O));
		polyPatterns["BAAAA"]+=(P1*(1.0-P2)*(1.0-P3a)*(1.0-P3b)*(1.0-O));
	}
}

double Dfoil::getDFO()
{
	return DFO;
}

double Dfoil::getDIL()
{
	return DIL;
}

double Dfoil::getDFI()
{
	return DFI;
}

double Dfoil::getDOL()
{
	return DOL;
}

double Dfoil::getChiSqrDFO()
{
	return xsqrDFO;
}

double Dfoil::getChiSqrDIL()
{
	return xsqrDIL;
}

double Dfoil::getChiSqrDFI()
{
	return xsqrDFI;
}

double Dfoil::getChiSqrDOL()
{
	return xsqrDOL;
}

int Dfoil::getPatternInt(std::string pattern)
{
	std::unordered_map<std::string, int>::const_iterator got=patterns.find(pattern);
	if(got!=patterns.end())
	{
		return patterns[pattern];
	}
	else
	{
		return 0;
	}
}

void Dfoil::calcChiSqr()
{
	xsqrDFO = this->chisqr(patterns["BABAA"], patterns["BBBAA"], patterns["ABABA"], patterns["AAABA"], patterns["BAABA"], patterns["BBABA"], patterns["ABBAA"], patterns["AABAA"]);
	xsqrDIL = this->chisqr(patterns["ABBAA"], patterns["BBBAA"], patterns["BAABA"], patterns["AAABA"], patterns["ABABA"], patterns["BBABA"], patterns["BABAA"], patterns["AABAA"]);
	xsqrDFI = this->chisqr(patterns["BABAA"], patterns["BABBA"], patterns["ABABA"], patterns["ABAAA"], patterns["ABBAA"], patterns["ABBBA"], patterns["BAABA"], patterns["BAAAA"]);
	xsqrDOL = this->chisqr(patterns["BAABA"], patterns["BABBA"], patterns["ABBAA"], patterns["ABAAA"], patterns["ABABA"], patterns["ABBBA"], patterns["BABAA"], patterns["BAAAA"]);
}

void Dfoil::polyCalcChiSqr()
{
	xsqrDFO = this->chisqr(polyPatterns["BABAA"], polyPatterns["BBBAA"], polyPatterns["ABABA"], polyPatterns["AAABA"], polyPatterns["BAABA"], polyPatterns["BBABA"], polyPatterns["ABBAA"], polyPatterns["AABAA"]);
	xsqrDIL = this->chisqr(polyPatterns["ABBAA"], polyPatterns["BBBAA"], polyPatterns["BAABA"], polyPatterns["AAABA"], polyPatterns["ABABA"], polyPatterns["BBABA"], polyPatterns["BABAA"], polyPatterns["AABAA"]);
	xsqrDFI = this->chisqr(polyPatterns["BABAA"], polyPatterns["BABBA"], polyPatterns["ABABA"], polyPatterns["ABAAA"], polyPatterns["ABBAA"], polyPatterns["ABBBA"], polyPatterns["BAABA"], polyPatterns["BAAAA"]);
	xsqrDOL = this->chisqr(polyPatterns["BAABA"], polyPatterns["BABBA"], polyPatterns["ABBAA"], polyPatterns["ABAAA"], polyPatterns["ABABA"], polyPatterns["ABBBA"], polyPatterns["BABAA"], polyPatterns["BAAAA"]);
}

void Dfoil::calcStats(unsigned int length)
{
	calcDFI();
	calcDFO();
	calcDIL();
	calcDOL();

	numloci = length;

	std::cout << "Found " << length << " biallelic loci" << std::endl;
	std::cout << "AAABA = " << patterns["AAABA"] << std::endl;
	std::cout << "AABAA = " << patterns["AABAA"] << std::endl;
	std::cout << "ABAAA = " << patterns["ABAAA"] << std::endl;
	std::cout << "ABABA = " << patterns["ABABA"] << std::endl;
	std::cout << "ABBAA = " << patterns["ABBAA"] << std::endl;
	std::cout << "ABBBA = " << patterns["ABBBA"] << std::endl;
	std::cout << "BAAAA = " << patterns["BAAAA"] << std::endl;
	std::cout << "BAABA = " << patterns["BAABA"] << std::endl;
	std::cout << "BABAA = " << patterns["BABAA"] << std::endl;
	std::cout << "BABBA = " << patterns["BABBA"] << std::endl;
	std::cout << "BBABA = " << patterns["BBABA"] << std::endl;
	std::cout << "BBBAA = " << patterns["BBBAA"] << std::endl;
	std::cout << "DFO = " << DFO << std::endl;
	std::cout << "DIL = " << DIL << std::endl;
	std::cout << "DFI = " << DFI << std::endl;
	std::cout << "DOL = " << DOL << std::endl << std::endl;
}

void Dfoil::polyCalcStats(unsigned int length)
{
	polyCalcDFI();
	polyCalcDFO();
	polyCalcDIL();
	polyCalcDOL();

	numloci = length;

	std::cout << "Found " << length << " biallelic loci" << std::endl;
	std::cout << "AAABA = " << polyPatterns["AAABA"] << std::endl;
	std::cout << "AABAA = " << polyPatterns["AABAA"] << std::endl;
	std::cout << "ABAAA = " << polyPatterns["ABAAA"] << std::endl;
	std::cout << "ABABA = " << polyPatterns["ABABA"] << std::endl;
	std::cout << "ABBAA = " << polyPatterns["ABBAA"] << std::endl;
	std::cout << "ABBBA = " << polyPatterns["ABBBA"] << std::endl;
	std::cout << "BAAAA = " << polyPatterns["BAAAA"] << std::endl;
	std::cout << "BAABA = " << polyPatterns["BAABA"] << std::endl;
	std::cout << "BABAA = " << polyPatterns["BABAA"] << std::endl;
	std::cout << "BABBA = " << polyPatterns["BABBA"] << std::endl;
	std::cout << "BBABA = " << polyPatterns["BBABA"] << std::endl;
	std::cout << "BBBAA = " << polyPatterns["BBBAA"] << std::endl;
	std::cout << "DFO = " << DFO << std::endl;
	std::cout << "DIL = " << DIL << std::endl;
	std::cout << "DFI = " << DFI << std::endl;
	std::cout << "DOL = " << DOL << std::endl << std::endl;
}

void Dfoil::polyCalcStats()
{
	polyCalcDFI();
	polyCalcDFO();
	polyCalcDIL();
	polyCalcDOL();
}

void Dfoil::calcStats()
{
	calcDFI();
	calcDFO();
	calcDIL();
	calcDOL();
}

void Dfoil::bootstrap(int mpiboot, int bootstrap, std::unordered_map <std::string,int> &indlist, std::default_random_engine &generator, int ntaxa, locusfile &current, std::vector<int> &keep, int my_rank, int i, int combs, std::string *indarray, std::string output, bool hetIgnore, bool hetInclude)
{
	double *allBootDFO;
	double *allBootDIL;
	double *allBootDFI;
	double *allBootDOL;
	allBootDFO = new double[bootstrap];
	allBootDIL = new double[bootstrap];
	allBootDFI = new double[bootstrap];            
	allBootDOL = new double[bootstrap];
      
	//do bootstrapping
	double *bootDFO;
	double *bootDIL;
	double *bootDFI;
	double *bootDOL;
	bootDFO = new double[bootstrap];
	bootDIL = new double[bootstrap];
	bootDFI = new double[bootstrap];
	bootDOL = new double[bootstrap];

	bootproc(mpiboot, keep, bootDFO, bootDIL, bootDFI, bootDOL, current, indlist, generator, ntaxa, hetIgnore, hetInclude, my_rank); //call private bootstrap procedure within this class

	MPI_Barrier(MPI_COMM_WORLD); //barrier after bootstrapping
	
	MPI_Gather(bootDFO, mpiboot, MPI_DOUBLE, allBootDFO, mpiboot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(bootDIL, mpiboot, MPI_DOUBLE, allBootDIL, mpiboot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(bootDFI, mpiboot, MPI_DOUBLE, allBootDFI, mpiboot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(bootDOL, mpiboot, MPI_DOUBLE, allBootDOL, mpiboot, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
		
	if(my_rank == 0)
	{
		std::cout << "statistics calculated for test " << i+1 << " of " << combs << std::endl;
		//calculate avg from bootstrapping
		avgDFO = this->average(allBootDFO, bootstrap);
		avgDIL = this->average(allBootDIL, bootstrap);
		avgDFI = this->average(allBootDFI, bootstrap);
		avgDOL = this->average(allBootDOL, bootstrap);

		//calculate standard deviation
		sdDFO = this->stdev(allBootDFO, avgDFO, bootstrap);
		sdDIL = this->stdev(allBootDIL, avgDIL, bootstrap);
		sdDFI = this->stdev(allBootDFI, avgDFI, bootstrap);
		sdDOL = this->stdev(allBootDOL, avgDOL, bootstrap);

		//calculate Z scores
		ZDFO = this->calcZ(DFO, sdDFO);
		ZDIL = this->calcZ(DIL, sdDIL);
		ZDFI = this->calcZ(DFI, sdDFI);
		ZDOL = this->calcZ(DOL, sdDOL);

		//calculate p-values for Z scores
		boost::math::normal_distribution<> zdist(0.0, 1.0);
		ZDFOpval = 2.0*(1-boost::math::cdf(zdist, abs(ZDFO)));
		ZDILpval = 2.0*(1-boost::math::cdf(zdist, abs(ZDIL)));
		ZDFIpval = 2.0*(1-boost::math::cdf(zdist, abs(ZDFI)));
		ZDOLpval = 2.0*(1-boost::math::cdf(zdist, abs(ZDOL)));

		//calculate chi squared values
		if(hetIgnore==true || (hetIgnore==false && hetInclude==false))
		{
			calcChiSqr();
		}
		else if(hetInclude==true)
		{
			polyCalcChiSqr();
		}
		else
		{
			std::cerr << "This code should not be reachable." << std::endl;
			exit(EXIT_FAILURE);
		}


		//calculate chi squared p-values

		int df = 1;
		boost::math::chi_squared Xdist(df);

		XpvalDFO = 1 - boost::math::cdf(Xdist, xsqrDFO);
		XpvalDIL = 1 - boost::math::cdf(Xdist, xsqrDIL);
		XpvalDFI = 1 - boost::math::cdf(Xdist, xsqrDFI);
		XpvalDOL = 1 - boost::math::cdf(Xdist, xsqrDOL);
			
		std::cout << "Chi-squared test calculated" << std::endl;
			
		//print out statistics
		writeout(indarray, output, i, hetIgnore, hetInclude );
		std::cout << "D statistics written" << std::endl;
		for(int i=0; i<5; i++)
		{
			std::cout << indarray[i] << std::endl;
		}
		std::cout << std::endl;
	}
	
	delete[] indarray;
	delete[] allBootDFO;
	delete[] allBootDIL;
	delete[] allBootDFI;
	delete[] allBootDOL;
	delete[] bootDFO;
	delete[] bootDIL;
	delete[] bootDFI;
	delete[] bootDOL;     
}

void Dfoil::bootproc(int bootstrap, std::vector<int> &keep, double *bootDFO, double *bootDIL, double *bootDFI, double *bootDOL, locusfile &file, std::unordered_map <std::string,int> &indlist, std::default_random_engine &generator, int ntaxa, bool hetIgnore, bool hetInclude, int my_rank)
{
	std::uniform_int_distribution<int> uniform(0, keep.size()-1);
	
	for(int i=0; i<bootstrap; i++)
	{
		std::vector<int> bootloci;
		//generate random ints
		for(unsigned int j=0; j<keep.size(); j++)
		{
			int temp = keep.at(uniform(generator));
			bootloci.push_back(temp);
		}
		
		fourtax bootrep(keep.size(), ntaxa);
		
		
		Dfoil bootDfoiltest;
		//make switch to do poly test or other
		
		if(hetIgnore==true || (hetIgnore==false && hetInclude==false))
		{
			bootrep.populateDtest(bootloci, file, indlist, generator, my_rank, ntaxa);
			bootDfoiltest.calcDs(bootrep, keep.size(), ntaxa, file, my_rank);
			bootDfoiltest.calcStats();
		}
		else if(hetInclude==true)
		{
			bootrep.populateDtest(bootloci, file, indlist, my_rank, ntaxa);
			bootDfoiltest.calcPolyDs(bootrep, keep.size());
			bootDfoiltest.polyCalcStats();
		}
		
		
		bootDFO[i] = bootDfoiltest.getDFO();
		bootDIL[i] = bootDfoiltest.getDIL();
		bootDFI[i] = bootDfoiltest.getDFI();
		bootDOL[i] = bootDfoiltest.getDOL();
	    
	}
}

void Dfoil::write(std::string *array, std::ofstream &outfile, int i, bool hetIgnore, bool hetInclude)
{
        if(i==0)
        {
            outfile << "O" << "\t" 
                    << "P4" << "\t"
                    << "P3" << "\t" 
                    << "P2" << "\t" 
                    << "P1" << "\t" 
                    << "AAABA" << "\t" 
                    << "AABAA" << "\t" 
                    << "ABAAA" << "\t"
                    << "ABABA" << "\t"
                    << "ABBAA" << "\t"
                    << "ABBBA" << "\t"
                    << "BAAAA" << "\t"
                    << "BAABA" << "\t"
                    << "BABAA" << "\t"
                    << "BABBA" << "\t"
                    << "BBABA" << "\t"
                    << "BBBAA" << "\t"
		    << "nloci" << "\t"
                    << "DFO" << "\t" 
                    << "DIL" << "\t" 
                    << "DFI" << "\t"     
                    << "DOL" << "\t"
                    << "SD_DFO" << "\t" 
                    << "SD_DIL" << "\t" 
                    << "SD_DFI" << "\t"  
                    << "SD_DOL" << "\t"
                    << "X^2_DFO" << "\t" 
                    << "X^2_DIL" << "\t" 
                    << "X^2_DFI" << "\t"                     
                    << "X^2_DOL" << "\t"
                    << "X^2_p_DFO" << "\t"
                    << "X^2_p_DIL" << "\t" 
                    << "X^2_p_DFI" << "\t"     
                    << "X^2_p_DOL" << "\t"
                    << "Z_DFO" << "\t" 
                    << "Z_DIL" << "\t"
                    << "Z_DFI" << "\t"  
                    << "Z_DOL" << "\t"
                    << "Z_p_DFO" << "\t"
                    << "Z_p_DIL" << "\t"
                    << "Z_p_DFI" << "\t"
                    << "Z_p_DOL" << std::endl;
        }
        
        for(int j=0; j<5; j++)
        {
            outfile << array[j] << "\t";
        }
        
        if(hetIgnore==true || (hetIgnore==false && hetInclude==false))
	{
		outfile << patterns["AAABA"] << "\t";
		outfile << patterns["AABAA"] << "\t";
		outfile << patterns["ABAAA"] << "\t";
		outfile << patterns["ABABA"] << "\t";
		outfile << patterns["ABBAA"] << "\t";
		outfile << patterns["ABBBA"] << "\t";
		outfile << patterns["BAAAA"] << "\t";
		outfile << patterns["BAABA"] << "\t";
		outfile << patterns["BABAA"] << "\t";
		outfile << patterns["BABBA"] << "\t";
		outfile << patterns["BBABA"] << "\t";
		outfile << patterns["BBBAA"] << "\t";    
	}
	else if(hetInclude==true)
	{
		outfile << polyPatterns["AAABA"] << "\t";
		outfile << polyPatterns["AABAA"] << "\t";
		outfile << polyPatterns["ABAAA"] << "\t";
		outfile << polyPatterns["ABABA"] << "\t";
		outfile << polyPatterns["ABBAA"] << "\t";
		outfile << polyPatterns["ABBBA"] << "\t";
		outfile << polyPatterns["BAAAA"] << "\t";
		outfile << polyPatterns["BAABA"] << "\t";
		outfile << polyPatterns["BABAA"] << "\t";
		outfile << polyPatterns["BABBA"] << "\t";
		outfile << polyPatterns["BBABA"] << "\t";
		outfile << polyPatterns["BBBAA"] << "\t";    
	}
	else
	{
		std::cerr << "This code should not be reachable." << std::endl;
		exit(EXIT_FAILURE);
	}
        
      	outfile << numloci << "\t";
        outfile << std::fixed << std::setprecision(4) << DFO << "\t";
        outfile << std::fixed << std::setprecision(4) << DIL << "\t";
        outfile << std::fixed << std::setprecision(4) << DFI << "\t";
        outfile << std::fixed << std::setprecision(4) << DOL << "\t";        
        outfile << std::fixed << std::setprecision(4) << sdDFO << "\t";
        outfile << std::fixed << std::setprecision(4) << sdDIL << "\t";
        outfile << std::fixed << std::setprecision(4) << sdDFI << "\t";
        outfile << std::fixed << std::setprecision(4) << sdDOL << "\t";        
        outfile << std::fixed << std::setprecision(4) << xsqrDFO << "\t";
        outfile << std::fixed << std::setprecision(4) << xsqrDIL << "\t";
        outfile << std::fixed << std::setprecision(4) << xsqrDFI << "\t";
        outfile << std::fixed << std::setprecision(4) << xsqrDOL << "\t";        
        outfile << std::fixed << std::setprecision(4) << XpvalDFO << "\t";
        outfile << std::fixed << std::setprecision(4) << XpvalDIL << "\t";
        outfile << std::fixed << std::setprecision(4) << XpvalDFI << "\t";
        outfile << std::fixed << std::setprecision(4) << XpvalDOL << "\t";        
        outfile << std::fixed << std::setprecision(4) << ZDFO << "\t";
        outfile << std::fixed << std::setprecision(4) << ZDIL << "\t";
        outfile << std::fixed << std::setprecision(4) << ZDFI << "\t";
        outfile << std::fixed << std::setprecision(4) << ZDOL << "\t";        
        outfile << std::fixed << std::setprecision(4) << ZDFOpval << "\t";
        outfile << std::fixed << std::setprecision(4) << ZDILpval << "\t";
        outfile << std::fixed << std::setprecision(4) << ZDFIpval << "\t";
        outfile << std::fixed << std::setprecision(4) << ZDOLpval << std::endl;        
}
/*
void Dfoil::popCalcZ()
{
  
}
*/
