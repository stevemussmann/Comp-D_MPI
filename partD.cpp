/* 
 * File:   partD.cpp
 * Author: Steve
 * 
 * Created on July 23, 2015, 8:22 PM
 */

#include "partD.h"
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

partD::partD() 
{
    D12=0;
    D1=0;
    D2=0;
    xsqrD1 = 0.0;
    xsqrD2 = 0.0;
    xsqrD12 = 0.0;
    patterns["ABBBA"] = 0;
    patterns["BABBA"] = 0;
    patterns["ABBAA"] = 0;
    patterns["BABAA"] = 0;
    patterns["ABABA"] = 0;
    patterns["BAABA"] = 0;
    
    //avg values from bootstrapping
    avgD1 = 0.0;
    avgD2 = 0.0;
    avgD12 = 0.0;
    
    //stdev values from bootstrapping
    sdD1 = 0.0;
    sdD2 = 0.0;
    sdD12 = 0.0;
    
    //Z scores from bootstrapping
    ZD1 = 0.0;
    ZD2 = 0.0;
    ZD12 = 0.0;
    
    //p-values for Z scores
    ZD1pval = 0.0;
    ZD2pval = 0.0;
    ZD12pval = 0.0;
    
    //pvals for chi-squared test
    XpvalD1 = 0.0;
    XpvalD2 = 0.0;
    XpvalD12 = 0.0;
}

partD::partD(const partD& orig) {
}

partD::~partD() {
}

void partD::calcD12()
{
    int denom=patterns["ABBBA"]+patterns["BABBA"];
    if(denom == 0)
    {
        D12=0.0;
    }
    else
    {    
        int num=patterns["ABBBA"]-patterns["BABBA"];
        D12=(double)num/(double)denom;
    }
}

void partD::calcD1()
{
    int denom = patterns["ABBAA"] + patterns["BABAA"];
    if(denom == 0)
    {
        D1=0.0;
    }
    else
    {
        int num = patterns["ABBAA"]-patterns["BABAA"];
        
        D1=(double)num/(double)denom;
    }
}

void partD::calcD2()
{
    int denom=patterns["ABABA"]+patterns["BAABA"];
    if(denom == 0)
    {
        D2=0.0;
    }
    else
    {    
        int num=patterns["ABABA"]-patterns["BAABA"];
        
        D2=(double)num/(double)denom;
    }
}

void partD::polyCalcD12()
{
    double denom=polyPatterns["ABBBA"]+polyPatterns["BABBA"];
    if( (int)denom == 0)
    {
        D12=0.0;
    }
    else
    {    
        double num=polyPatterns["ABBBA"]-polyPatterns["BABBA"];
        D12=num/denom;
    }
}

void partD::polyCalcD1()
{
    double denom = polyPatterns["ABBAA"] + polyPatterns["BABAA"];
    if( (int)denom == 0)
    {
        D1=0.0;
    }
    else
    {
        double num = polyPatterns["ABBAA"]-polyPatterns["BABAA"];
        
        D1=num/denom;
    }
}

void partD::polyCalcD2()
{
    double denom=polyPatterns["ABABA"]+polyPatterns["BAABA"];
    if( (int)denom == 0)
    {
        D2=0.0;
    }
    else
    {    
        double num=polyPatterns["ABABA"]-polyPatterns["BAABA"];
        
        D2=num/denom;
    }
}

void partD::calcDs(fourtax &dtest, unsigned int length, int ntaxa, locusfile &file, int my_rank)
{
    for(unsigned int i=0; i<length; i++)
    {
        dtest.calculatePattern(i, ntaxa, file, my_rank);
        patterns[dtest.getPattern(i)]++;
    }
}

void partD::calcPolyDs(fourtax &dtest, int nloci)
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
	}

}

double partD::getD1()
{
    return D1;
}

double partD::getD2()
{
    return D2;
}

double partD::getD12()
{
    return D12;
}

double partD::getChiSqrD1()
{
    return xsqrD1;
}

double partD::getChiSqrD2()
{
    return xsqrD2;
}

double partD::getChiSqrD12()
{
    return xsqrD12;
}

int partD::getPatternInt(std::string pattern)
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

void partD::calcChiSqr()
{
    xsqrD1 = this->chisqr(patterns["ABBAA"], patterns["BABAA"]);
    xsqrD2 = this->chisqr(patterns["ABABA"], patterns["BAABA"]);
    xsqrD12 = this->chisqr(patterns["ABBBA"], patterns["BABBA"]);
}

void partD::polyCalcChiSqr()
{
    xsqrD1 = this->chisqr(polyPatterns["ABBAA"], polyPatterns["BABAA"]);
    xsqrD2 = this->chisqr(polyPatterns["ABABA"], polyPatterns["BAABA"]);
    xsqrD12 = this->chisqr(polyPatterns["ABBBA"], polyPatterns["BABBA"]);
}

void partD::calcStats(unsigned int length)
{
    calcD1();
    calcD2();
    calcD12();
    
    std::cout << "Found " << length << " biallelic loci" << std::endl;
    std::cout << "ABBBA = " << patterns["ABBBA"] << std::endl;
    std::cout << "BABBA = " << patterns["BABBA"] << std::endl;
    std::cout << "ABBAA = " << patterns["ABBAA"] << std::endl;
    std::cout << "BABAA = " << patterns["BABAA"] << std::endl;
    std::cout << "ABABA = " << patterns["ABABA"] << std::endl;
    std::cout << "BAABA = " << patterns["BAABA"] << std::endl;
    std::cout << "D1 = " << D1 << std::endl;
    std::cout << "D2 = " << D2 << std::endl;
    std::cout << "D12 = " << D12 << std::endl << std::endl;
}

void partD::polyCalcStats(unsigned int length)
{
    polyCalcD1();
    polyCalcD2();
    polyCalcD12();
    
    std::cout << "Found " << length << " biallelic loci" << std::endl;
    std::cout << "ABBBA = " << polyPatterns["ABBBA"] << std::endl;
    std::cout << "BABBA = " << polyPatterns["BABBA"] << std::endl;
    std::cout << "ABBAA = " << polyPatterns["ABBAA"] << std::endl;
    std::cout << "BABAA = " << polyPatterns["BABAA"] << std::endl;
    std::cout << "ABABA = " << polyPatterns["ABABA"] << std::endl;
    std::cout << "BAABA = " << polyPatterns["BAABA"] << std::endl;
    std::cout << "D1 = " << D1 << std::endl;
    std::cout << "D2 = " << D2 << std::endl;
    std::cout << "D12 = " << D12 << std::endl << std::endl;
}

void partD::polyCalcStats()
{
    polyCalcD1();
    polyCalcD2();
    polyCalcD12();
}

void partD::calcStats()
{
    calcD1();
    calcD2();
    calcD12();
}

void partD::bootstrap(int mpiboot, int bootstrap, std::unordered_map <std::string,int> &indlist, 
		      std::default_random_engine &generator, int ntaxa, locusfile &current, 
		      std::vector<int> &keep, int my_rank, int i, int combs, std::string *indarray, 
		      std::string output, bool hetIgnore, bool hetInclude)
{
    double *allBootD1;

    double *allBootD2;

    double *allBootD12;

    allBootD1 = new double[bootstrap];

    allBootD2 = new double[bootstrap];

    allBootD12 = new double[bootstrap];
  
    //do bootstrapping
    double *bootD1;

    double *bootD2;

    double *bootD12;

    bootD1 = new double[mpiboot];

    bootD2 = new double[mpiboot];

    bootD12 = new double[mpiboot];

    bootproc(mpiboot, keep, bootD1, bootD2, bootD12, current, indlist, generator, ntaxa, hetIgnore, hetInclude, my_rank); //call private bootstrap procedure within this class

    MPI_Barrier(MPI_COMM_WORLD); //barrier after bootstrapping
    
    MPI_Gather(bootD1, mpiboot, MPI_DOUBLE, allBootD1, mpiboot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(bootD2, mpiboot, MPI_DOUBLE, allBootD2, mpiboot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(bootD12, mpiboot, MPI_DOUBLE, allBootD12, mpiboot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            
    if(my_rank == 0)
    {
	std::cout << "statistics calculated for test " << i+1 << " of " << combs << std::endl;
	//calculate avg from bootstrapping
        avgD1 = this->average(allBootD1, bootstrap);
        avgD2 = this->average(allBootD2, bootstrap);
        avgD12 = this->average(allBootD12, bootstrap);

        //calculate standard deviation

        sdD1 = this->stdev(allBootD1, avgD1, bootstrap);
        sdD2 = this->stdev(allBootD2, avgD2, bootstrap);
        sdD12 = this->stdev(allBootD12, avgD12, bootstrap);

        //calculate Z scores

        ZD1 = this->calcZ(D1, sdD1);
        ZD2 = this->calcZ(D2, sdD2);
        ZD12 = this->calcZ(D12, sdD12);

        //calculate p-values for Z scores

        boost::math::normal_distribution<> zdist(0.0, 1.0);
        ZD1pval = 2.0*(1-boost::math::cdf(zdist, abs(ZD1)));
        ZD2pval = 2.0*(1-boost::math::cdf(zdist, abs(ZD2)));
        ZD12pval = 2.0*(1-boost::math::cdf(zdist, abs(ZD12)));

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

        XpvalD1 = 1 - boost::math::cdf(Xdist, xsqrD1);

        XpvalD2 = 1 - boost::math::cdf(Xdist, xsqrD2);

        XpvalD12 = 1 - boost::math::cdf(Xdist, xsqrD12);
                
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

    delete[] bootD1;

    delete[] bootD2;

    delete[] bootD12;

    delete[] allBootD1;

    delete[] allBootD2;

    delete[] allBootD12;
}

void partD::bootproc(int bootstrap, std::vector<int> &keep, double *bootD1, double *bootD2, double *bootD12, locusfile &file, std::unordered_map <std::string,int> &indlist, std::default_random_engine &generator, int ntaxa, bool hetIgnore, bool hetInclude, int my_rank)
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
        partD bootpartDtest;
	
	if(hetIgnore==true || (hetIgnore==false && hetInclude==false))
	{
		bootrep.populateDtest(bootloci, file, indlist, generator, my_rank, ntaxa);
		bootpartDtest.calcDs(bootrep, keep.size(), ntaxa, file, my_rank);
		bootpartDtest.calcStats();
	}
	else if(hetInclude==true)
	{
		bootrep.populateDtest(bootloci, file, indlist, my_rank, ntaxa);
		bootpartDtest.calcPolyDs(bootrep, keep.size());
		bootpartDtest.polyCalcStats();
	}
	
	
        bootD1[i] = bootpartDtest.getD1();
        bootD2[i] = bootpartDtest.getD2();
        bootD12[i] = bootpartDtest.getD12();
	
    }
}

void partD::write(std::string *array, std::ofstream &outfile, int i, bool hetIgnore, bool hetInclude)
{
        if(i==0)
        {
	    outfile << "O" << "\t"
            << "P4" << "\t"
            << "P3" << "\t"
            << "P2" << "\t"
            << "P1" << "\t"
            << "ABBBA" << "\t"
            << "BABBA" << "\t"
            << "ABBAA" << "\t"
            << "BABAA" << "\t"
            << "ABABA" << "\t"
            << "BAABA" << "\t"
            << "D1" << "\t"
            << "D2" << "\t"
            << "D12" << "\t"
            << "SD_D1" << "\t"
            << "SD_D2" << "\t"
            << "SD_D12" << "\t"
            << "X^2_D1" << "\t"
            << "X^2_D2" << "\t"
            << "X^2_D12" << "\t"
            << "X^2_p_D1" << "\t"
            << "X^2_p_D2" << "\t"
            << "X^2_p_D12" << "\t"
            << "Z_D1" << "\t"
            << "Z_D2" << "\t"
            << "Z_D12" << "\t"
            << "Z_p_D1" << "\t"
            << "Z_p_D2" << "\t"
            << "Z_p_D12" << std::endl;
        }
        
        for(int j=0; j<5; j++)
        {
            outfile << array[j] << "\t";
        }
        
        if(hetIgnore==true || (hetIgnore==false && hetInclude==false))
	{
		outfile << patterns["ABBBA"] << "\t";
		outfile << patterns["BABBA"] << "\t";
		outfile << patterns["ABBAA"] << "\t";
		outfile << patterns["BABAA"] << "\t";
		outfile << patterns["ABABA"] << "\t";
		outfile << patterns["BAABA"] << "\t";
	}
	else if(hetInclude==true)
	{
		outfile << polyPatterns["ABBBA"] << "\t";
		outfile << polyPatterns["BABBA"] << "\t";
		outfile << polyPatterns["ABBAA"] << "\t";
		outfile << polyPatterns["BABAA"] << "\t";
		outfile << polyPatterns["ABABA"] << "\t";
		outfile << polyPatterns["BAABA"] << "\t";
	}
	else
	{
		std::cerr << "This code should not be reachable." << std::endl;
		exit(EXIT_FAILURE);
	}

        outfile << std::fixed << std::setprecision(4) << D1 << "\t";
        outfile << std::fixed << std::setprecision(4) << D2 << "\t";
        outfile << std::fixed << std::setprecision(4) << D12 << "\t";
        outfile << std::fixed << std::setprecision(4) << sdD1 << "\t";
        outfile << std::fixed << std::setprecision(4) << sdD2 << "\t";
        outfile << std::fixed << std::setprecision(4) << sdD12 << "\t";
        outfile << std::fixed << std::setprecision(4) << xsqrD1 << "\t";
        outfile << std::fixed << std::setprecision(4) << xsqrD2 << "\t";
        outfile << std::fixed << std::setprecision(4) << xsqrD12 << "\t";
        outfile << std::fixed << std::setprecision(4) << XpvalD1 << "\t";
        outfile << std::fixed << std::setprecision(4) << XpvalD2 << "\t";
        outfile << std::fixed << std::setprecision(4) << XpvalD12 << "\t";
        outfile << std::fixed << std::setprecision(4) << ZD1 << "\t";
        outfile << std::fixed << std::setprecision(4) << ZD2 << "\t";
        outfile << std::fixed << std::setprecision(4) << ZD12 << "\t";
        outfile << std::fixed << std::setprecision(4) << ZD1pval << "\t";
        outfile << std::fixed << std::setprecision(4) << ZD2pval << "\t";
        outfile << std::fixed << std::setprecision(4) << ZD12pval << std::endl;
}
/*
void partD::popCalcZ()
{
  
}
*/