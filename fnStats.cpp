#include "fnFiles.h"
#include "fnStats.h"

//#include <algorithm>
//#include <cstdlib>
//#include <fstream>
//#include <iostream>
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
}

void fnStats::calcAllFreqs()
{
	calcFreqs(Afreqs);
	calcFreqs(Bfreqs);
	calcFreqs(Cfreqs);
	calcFreqs(Dfreqs);
}

void fnStats::calcFreqs(std::vector<std::unordered_map <std::string,int> > &vec)
{

}

