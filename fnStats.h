#ifndef FNSTATS_H
#define FNSTATS_H

//#include "fnFiles.h"

#include <string>
#include <unordered_map>
#include <vector>

class fnStats {
	public:
		fnStats(int l);
		void calcAllFreqs();
	private:
		//record alleles as anc (ancestral) and derived (der)
		std::vector<std::unordered_map <std::string,int> > Afreqs;
		std::vector<std::unordered_map <std::string,int> > Bfreqs;
		std::vector<std::unordered_map <std::string,int> > Cfreqs;
		std::vector<std::unordered_map <std::string,int> > Dfreqs;
		std::unordered_map<int,int> blacklist;
		void calcFreqs(std::vector<std::unordered_map <std::string,int> > &vec);
};

#endif
