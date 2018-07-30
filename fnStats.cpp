#include "fnFiles.h"
#include "fnStats.h"

//#include <algorithm>
//#include <cstdlib>
//#include <fstream>
#include <iostream>
#include <iterator>
//#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

fnStats::fnStats(int l)
{
	Afreqs.resize(l);
	Bfreqs.resize(l);
	Cfreqs.resize(l);
	Dfreqs.resize(l);
	ancestral.resize(l);
}

void fnStats::calcAllFreqs()
{
	calcFreqs(Afreqs);
	calcFreqs(Bfreqs);
	calcFreqs(Cfreqs);
	calcFreqs(Dfreqs);
}

void fnStats::findAncestral(fnFiles &f)
{
	for(unsigned int i=0; i<ancestral.size(); i++)
	{
		std::unordered_map<std::string,int> l = f.getOutgroupLocus(i); //get locus for outgroup

		std::unordered_map<std::string,int>::iterator it = l.begin(); //start iterator
		if(l.size() == 1)
		{
			ancestral[i] = it->first; //if only one allele, set it as ancestral state.
		}
		else
		{
			std::string anc;
			int ancCount = 0;
			while(it != l.end())
			{
				if(it->second > ancCount){
					ancCount = it->second;
					anc = it->first;
				}
				it++;
			}
			ancestral[i] = anc;
		}
		
		std::cout << ancestral[i] << std::endl;
		//std::cout << l.size() << std::endl;
	}
}

void fnStats::calcFreqs(std::vector<std::unordered_map <std::string,int> > &vec)
{

}

