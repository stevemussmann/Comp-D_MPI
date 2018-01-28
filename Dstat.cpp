/* 
 * File:   Dstat.cpp
 * Author: Steve
 * 
 * Created on July 26, 2015, 11:48 AM
 */

#include "Dstat.h"
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

Dstat::Dstat() 
{
    D=0.0;
    ABBA=0;
    BABA=0;
    polyABBA=0.0;
    polyBABA=0.0;
    chisq=0.0;
    avg=0.0;
    sd=0.0;
    Z=0.0;
    Zpval=0.0;
    Xpval=0.0;
    numloci=0;
}

void Dstat::calcStats(unsigned int length)
{
    if(ABBA + BABA == 0)
    {
        D=0.0;
    }
    else
    {
        D = (double)(ABBA - BABA) / (double)(ABBA + BABA);
    }

    numloci = length;
    
    std::cout << "Found " << length << " biallelic loci" << std::endl;
    std::cout << "ABBA = " << ABBA << std::endl;
    std::cout << "BABA = " << BABA << std::endl << std::endl;
}

void Dstat::polyCalcStats(unsigned int length)
{
	if((int)(polyABBA+polyBABA)==0)
	{
		D=0.0;
	}
	else
	{
		D = (polyABBA - polyBABA) / (polyABBA + polyBABA);
	}

	numloci = length;
	
	std::cout << "Found " << length << " biallelic loci" << std::endl;
	std::cout << "ABBA = " << polyABBA << std::endl;
	std::cout << "BABA = " << polyBABA << std::endl << std::endl;
}

void Dstat::polyCalcStats()
{
	if((int)(polyABBA+polyBABA)==0)
	{
		D=0.0;
	}
	else
	{
		D = (polyABBA - polyBABA) / (polyABBA + polyBABA);
	}
}

void Dstat::calcStats()
{
    if(ABBA + BABA == 0)
    {
        D=0.0;
    }
    else
    {
        D = (double)(ABBA - BABA) / (double)(ABBA + BABA);
    }
}

void Dstat::calcDs(fourtax &dtest, unsigned int length, int ntaxa, locusfile &file, int my_rank)
{
	for(unsigned int i=0; i<length; i++)
	{
		//check if all taxa are present
		dtest.calculatePattern(i, ntaxa, file, my_rank);
		if(dtest.getPattern(i) == "ABBA")
		{
			ABBA++;
		}
		else if(dtest.getPattern(i) == "BABA")
		{
			BABA++;
		}
	}
}

void Dstat::calcPolyDs(fourtax &dtest, int nloci)
{
  	for(int i=0; i<nloci; i++)
	{
		double P1=dtest.getFreq(i, 3);
		double P2=dtest.getFreq(i, 2);
		double P3=dtest.getFreq(i, 1);
		double O=dtest.getFreq(i, 0);
	
		polyABBA+=((1.0-P1)*P2*P3*(1.0-O));
		polyBABA+=(P1*(1.0-P2)*P3*(1.0-O));
	}
}

double Dstat::getD()
{
    return D;
}

int Dstat::getABBA()
{
    return ABBA;
}

int Dstat::getBABA()
{
    return BABA;
}

double Dstat::getChiSqr()
{
    return chisq;
}

Dstat::Dstat(const Dstat& orig) {
}

Dstat::~Dstat() {
}

void Dstat::calcChiSqr() 
{
    chisq = this->chisqr(ABBA, BABA);
}

void Dstat::polyCalcChiSqr() 
{
    chisq = this->chisqr(polyABBA, polyBABA);
}

void Dstat::bootstrap(int mpiboot, int bootstrap, std::unordered_map <std::string,int> &indlist, 
		std::default_random_engine &generator, int ntaxa, locusfile &current, 
		std::vector<int> &keep, int my_rank, int i, int combs, std::string *indarray, 
		std::string output, bool hetIgnore, bool hetInclude)
{
    double *allBoot;
    allBoot = new double[bootstrap];
    //do bootstrapping
    double *bootD;
    bootD = new double[mpiboot];
    bootproc(mpiboot, keep, bootD, current, indlist, generator, ntaxa, hetIgnore, hetInclude, my_rank);

    MPI_Barrier(MPI_COMM_WORLD); //barrier after completing bootstrapping

    MPI_Gather(bootD, mpiboot, MPI_DOUBLE, allBoot, mpiboot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            
    if(my_rank == 0)
    {
	std::cout << "statistics calculated for test " << i+1 << " of " << combs << std::endl;
        //calculate avg and stdev from bootstrapping
        avg = this->average(allBoot, bootstrap);
        sd = this->stdev(allBoot, avg, bootstrap);
    
        //calculate Z score and pval
        Z = this->calcZ(D, sd);

        boost::math::normal_distribution<> zdist(0.0, 1.0);
        Zpval = 2.0*(1-boost::math::cdf(zdist, abs(Z)));
    
        //calculate chi squared and pval
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
        int df = 1;
        boost::math::chi_squared Xdist(df);
        Xpval = 1 - boost::math::cdf(Xdist, chisq);
                
        std::cout << "Chi-squared test calculated" << std::endl;
                
        //print out statistics
        writeout(indarray, output, i, hetIgnore, hetInclude );
        std::cout << "D statistics written" << std::endl;
        for(int i=0; i<4; i++)
        {
	    std::cout << indarray[i] << std::endl;
        }
        std::cout << std::endl;
    }
    delete[] indarray;
    delete[] bootD;
    delete[] allBoot;
}

void Dstat::bootproc(int bootstrap, std::vector<int> &keep, double *bootD, locusfile &file, 
		     std::unordered_map <std::string,int> &indlist, std::default_random_engine &generator, 
		     int ntaxa, bool hetIgnore, bool hetInclude, int my_rank)
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
        bootrep.populateDtest(bootloci, file, indlist, generator, my_rank, ntaxa);
        
        Dstat bootDtest;
        bootDtest.calcDs(bootrep, keep.size(), ntaxa, file, my_rank);
        bootDtest.calcStats();
        
        bootD[i] = bootDtest.getD();
    }
}

void Dstat::write(std::string *array, std::ofstream &outfile, int i, bool hetIgnore, bool hetInclude)
{
        if(i==0)
        {
            outfile << "O" << "\t" << "P3" << "\t" << "P2" << "\t" << "P1" << "\t" << "ABBA" << 
                "\t" << "BABA" << "\t" "nloci" << "\t" << "D" << "\t" << "STDEV" << "\t" << "X^2" << "\t" <<
                "X^2_pval" << "\t" << "Z-score" << "\t" << "Z pval" << std::endl;
        }
        
        for(int i=0; i<4; i++)
        {
            outfile << array[i] << "\t";
        }
        
        if(hetIgnore==true || (hetIgnore==false && hetInclude==false))
	{
		outfile << ABBA << "\t";
		outfile << BABA << "\t";
	}
	else if(hetInclude==true)
	{
		outfile << polyABBA << "\t";
		outfile << polyBABA << "\t";
	}
	else
	{
		std::cerr << "This code should not be reachable." << std::endl;
		exit(EXIT_FAILURE);
	}
        
	outfile << numloci << "\t";
        outfile << std::fixed << std::setprecision(4) << D << "\t";
        outfile << std::fixed << std::setprecision(4) << sd << "\t";
        outfile << std::fixed << std::setprecision(4) << chisq << "\t";
        outfile << std::fixed << std::setprecision(4) << Xpval << "\t";
        outfile << std::fixed << std::setprecision(4) << Z << "\t";
        outfile << std::fixed << std::setprecision(4) << Zpval << std::endl;
}
/*
void Dstat::popCalcZ()
{
  
}
*/
