#ifndef FNSTATS_H
#define FNSTATS_H

//#include "fnFiles.h"

#include "fnFiles.h"

#include <string>
#include <unordered_map>
#include <vector>

class fnStats {
	public:
		fnStats(int l, std::unordered_map<std::string,std::string> m);
		void calcAllFreqs(fnFiles &f, std::unordered_map<std::string,std::string> m);
		void findAncestral(fnFiles &f);
		void calcFstats(fnFiles &f);
	private:
		//record alleles as anc (ancestral) and derived (der)
		std::unordered_map<std::string,std::vector<std::unordered_map<std::string,double> > > freqs;
		std::vector<std::string> ancestral;
		std::unordered_map<int,int> blacklist;
		double hz(std::vector<int> &v, int t);
		double f2(double freqa, double freqb, double ha, double hb, int tota, int totb);
		double f3(double freqa, double freqb, double freqc, double hc, int totc);
		double f4(double freqa, double freqb, double freqc, double freqd);
		double fst(double freqa, double freqb);
};

#endif
