#include "popZpartD.h"
#include "partD.h"
#include "DstatParent.h"

#include "popZParent.h"

#include <iostream>

#include <boost/math/distributions/normal.hpp>

popZpartD::popZpartD(){

}

void popZpartD::add(DstatParent *d){
    partD* nd = dynamic_cast<partD*>(d);
    D1.push_back(nd->getD1());
    D2.push_back(nd->getD2());
    D12.push_back(nd->getD12());
}

void popZpartD::calcStats(std::string filename){
    partdZ(filename);
}

void popZpartD::partdZ(std::string filename){
	//calculate avg from bootstrapping
	double* D1arr = this->toArr(D1, D1.size());
	double* D2arr = this->toArr(D2, D2.size());
	double* D12arr = this->toArr(D12, D12.size());
	
	double avgD1 = this->average(D1arr, D1.size());
	double avgD2 = this->average(D2arr, D2.size());
	double avgD12 = this->average(D12arr, D12.size());
	
	//calculate standard deviation
	double sdD1 = this->stdev(D1arr, avgD1, D1.size());
	double sdD2 = this->stdev(D2arr, avgD2, D2.size());
	double sdD12 = this->stdev(D12arr, avgD12, D12.size());
	
	//calculate Z scores
	double ZD1 = this->calcZ(avgD1, sdD1);
	double ZD2 = this->calcZ(avgD2, sdD2);
	double ZD12 = this->calcZ(avgD12, sdD12);
	
	//calculate p-values for Z scores
	boost::math::normal_distribution<> zdist(0.0, 1.0);
	double ZD1pval = 2.0*(1-boost::math::cdf(zdist, fabs(ZD1)));
	double ZD2pval = 2.0*(1-boost::math::cdf(zdist, fabs(ZD2)));
	double ZD12pval = 2.0*(1-boost::math::cdf(zdist, fabs(ZD12)));

	std::ofstream popout;
	popout.open(filename, std::ios::out);
	if(popout.is_open())
	{
		popout << "Statistic" << "\t" << "Z-score" << "\t" << "P-val" << std::endl;
		popout << "D1" << "\t" << ZD1 << "\t" << ZD1pval << std::endl;
		popout << "D2" << "\t" << ZD2 << "\t" << ZD2pval << std::endl;
		popout << "D12" << "\t" << ZD12 << "\t" << ZD12pval << std::endl;
	}
	popout.close();
}
