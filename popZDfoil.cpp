#include "popZDfoil.h"
#include "Dfoil.h"
#include "DstatParent.h"

#include "popZParent.h"

#include <iostream>

#include <boost/math/distributions/normal.hpp>

popZDfoil::popZDfoil(){

}

void popZDfoil::add(DstatParent *d){
	Dfoil* nd = dynamic_cast<Dfoil*>(d);
	DFO.push_back(nd->getDFO());
	DIL.push_back(nd->getDIL());
	DFI.push_back(nd->getDFI());
	DOL.push_back(nd->getDOL());
}

void popZDfoil::calcStats(std::string filename){
	dfoilZ(filename);
}

void popZDfoil::dfoilZ(std::string filename){
	//calculate avg from bootstrapping
	double* DFOarr = this->toArr(DFO, DFO.size());
	double* DILarr = this->toArr(DIL, DIL.size());
	double* DFIarr = this->toArr(DFI, DFI.size());
	double* DOLarr = this->toArr(DOL, DOL.size());

	double avgDFO = this->average(DFOarr, DFO.size());
	double avgDIL = this->average(DILarr, DIL.size());
	double avgDFI = this->average(DFIarr, DFI.size());
	double avgDOL = this->average(DOLarr, DOL.size());

	//calculate standard deviation
	double sdDFO = this->stdev(DFOarr, avgDFO, DFO.size());
	double sdDIL = this->stdev(DILarr, avgDIL, DIL.size());
	double sdDFI = this->stdev(DFIarr, avgDFI, DFI.size());
	double sdDOL = this->stdev(DOLarr, avgDOL, DOL.size());

	//calculate Z scores
	double ZDFO = this->calcZ(avgDFO, sdDFO);
	double ZDIL = this->calcZ(avgDIL, sdDIL);
	double ZDFI = this->calcZ(avgDFI, sdDFI);
	double ZDOL = this->calcZ(avgDOL, sdDOL);

	//calculate p-values for Z scores
	boost::math::normal_distribution<> zdist(0.0, 1.0);
	double ZDFOpval = 2.0*(1-boost::math::cdf(zdist, abs(ZDFO)));
	double ZDILpval = 2.0*(1-boost::math::cdf(zdist, abs(ZDIL)));
	double ZDFIpval = 2.0*(1-boost::math::cdf(zdist, abs(ZDFI)));
	double ZDOLpval = 2.0*(1-boost::math::cdf(zdist, abs(ZDOL)));
   
	std::ofstream popout;
	popout.open(filename, std::ios::out);
	if(popout.is_open())
	{
		popout << "Statistic" << "\t" << "Z-score"  << "\t" << "P-val" << std::endl;
		popout << "DFO" << "\t" << ZDFO << "\t" << ZDFOpval << std::endl;
		popout << "DIL" << "\t" << ZDIL << "\t" << ZDILpval << std::endl;
		popout << "DFI" << "\t" << ZDFI << "\t" << ZDFIpval << std::endl;
		popout << "DOL" << "\t" << ZDOL << "\t" << ZDOLpval << std::endl;
	}
	popout.close();
}

