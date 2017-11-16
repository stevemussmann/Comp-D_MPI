#include "popZDstat.h"
#include "Dstat.h"
#include "DstatParent.h"

#include "popZParent.h"

#include <iostream>

#include <boost/math/distributions/normal.hpp>

popZDstat::popZDstat(){

}

void popZDstat::add(DstatParent *d){
    Dstat* nd = dynamic_cast<Dstat*>(d);
    D.push_back(nd->getD());
}

void popZDstat::calcStats(){
    dstatZ();
}

void popZDstat::dstatZ(){
	double* Darr = this->toArr(D, D.size());
	
	double avgD = this->average(Darr, D.size());
	
	double sdD = this->stdev(Darr, avgD, D.size());
	
	double ZD = this->calcZ(avgD, sdD);
	
	boost::math::normal_distribution<> zdist(0.0, 1.0);
	double ZDpval = 2.0*(1-boost::math::cdf(zdist, abs(ZD)));
	
	std::cout << ZD << "\t" << ZDpval << std::endl;
}