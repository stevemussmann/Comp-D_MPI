#include "fnFiles.h"
#include "fnStats.h"

//#include <algorithm>
//#include <cstdlib>
//#include <fstream>
#include <iostream>
//#include <iterator>
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
	for(int i=0; i<ancestral.size(); i++)
	{
		std::unordered_map <std::string,int> l = f.getOutgroupLocus(i);
		std::cout << l.size() << std::endl;
	}
}

void fnStats::calcFreqs(std::vector<std::unordered_map <std::string,int> > &vec)
{

}

