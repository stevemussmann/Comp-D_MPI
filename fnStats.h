#ifndef FNSTATS_H
#define FNSTATS_H

//#include "fnFiles.h"

#include "fnFiles.h"

#include <string>
#include <unordered_map>
#include <vector>

class fnStats {
	public:
		fnStats(int l);
		void calcAllFreqs(fnFiles &f);
		void findAncestral(fnFiles &f);
	private:
		//record alleles as anc (ancestral) and derived (der)
		std::vector<std::unordered_map <std::string,double> > Afreqs;
		std::vector<std::unordered_map <std::string,double> > Bfreqs;
		std::vector<std::unordered_map <std::string,double> > Cfreqs;
		std::vector<std::unordered_map <std::string,double> > Dfreqs;
		std::vector<std::string> ancestral;
		std::unordered_map<int,int> blacklist;
		void calcAFreqs(fnFiles &f);
		void calcBFreqs(fnFiles &f);
		void calcCFreqs(fnFiles &f);
		void calcDFreqs(fnFiles &f);
};

#endif
