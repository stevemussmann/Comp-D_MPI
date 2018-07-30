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

void fnStats::calcAllFreqs(fnFiles &f)
{
	calcAFreqs(f);
	calcBFreqs(f);
	calcCFreqs(f);
	calcDFreqs(f);
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
		
		//std::cout << ancestral[i] << std::endl;
		//std::cout << l.size() << std::endl;
	}
}

void fnStats::calcAFreqs(fnFiles &f)
{
	for(unsigned int i=0; i<Afreqs.size(); i++)
	{
		int total = 0;
		std::unordered_map<std::string,int> l = f.getALocus(i);
		std::unordered_map<std::string,int>::iterator it = l.begin();
		while(it != l.end())
		{
			total += it->second;
			it++;
		}

		it = l.begin();
		while(it != l.end())
		{
			double freq = it->second/(double)total;
			if(l.size() == 1){
				if(it->first == ancestral[i])
				{
					Afreqs[i]["anc"] = freq;
					Afreqs[i]["der"] = 0.0;
				}
				else if(it->first != ancestral[i])
				{
					Afreqs[i]["anc"] = 0.0;
					Afreqs[i]["der"] = freq;
				}
			}
			else if(l.size() == 2)
			{
				if(it->first == ancestral[i])
				{
					Afreqs[i]["anc"] = freq;
				}
				else
				{
					Afreqs[i]["der"] = freq;
				}
			}
			it++;
		}
		std::cout << i << std::endl;
		std::cout << "anc = " << Afreqs[i]["anc"] << std::endl;
		std::cout << "der = " << Afreqs[i]["der"] << std::endl;
	}
}

void fnStats::calcBFreqs(fnFiles &f)
{
	for(unsigned int i=0; i<Bfreqs.size(); i++)
	{
		int total = 0;
		std::unordered_map<std::string,int> l = f.getBLocus(i);
		std::unordered_map<std::string,int>::iterator it = l.begin();
		while(it != l.end())
		{
			total += it->second;
			it++;
		}

		it = l.begin();
		while(it != l.end())
		{
			double freq = it->second/(double)total;
			if(l.size() == 1){
				if(it->first == ancestral[i])
				{
					Bfreqs[i]["anc"] = freq;
					Bfreqs[i]["der"] = 0.0;
				}
				else if(it->first != ancestral[i])
				{
					Bfreqs[i]["anc"] = 0.0;
					Bfreqs[i]["der"] = freq;
				}
			}
			else if(l.size() == 2)
			{
				if(it->first == ancestral[i])
				{
					Bfreqs[i]["anc"] = freq;
				}
				else
				{
					Bfreqs[i]["der"] = freq;
				}
			}
			it++;
		}
		std::cout << i << std::endl;
		std::cout << "anc = " << Bfreqs[i]["anc"] << std::endl;
		std::cout << "der = " << Bfreqs[i]["der"] << std::endl;
	}
}

void fnStats::calcCFreqs(fnFiles &f)
{
	for(unsigned int i=0; i<Cfreqs.size(); i++)
	{
		int total = 0;
		std::unordered_map<std::string,int> l = f.getCLocus(i);
		std::unordered_map<std::string,int>::iterator it = l.begin();
		while(it != l.end())
		{
			total += it->second;
			it++;
		}

		it = l.begin();
		while(it != l.end())
		{
			double freq = it->second/(double)total;
			if(l.size() == 1){
				if(it->first == ancestral[i])
				{
					Cfreqs[i]["anc"] = freq;
					Cfreqs[i]["der"] = 0.0;
				}
				else if(it->first != ancestral[i])
				{
					Cfreqs[i]["anc"] = 0.0;
					Cfreqs[i]["der"] = freq;
				}
			}
			else if(l.size() == 2)
			{
				if(it->first == ancestral[i])
				{
					Cfreqs[i]["anc"] = freq;
				}
				else
				{
					Cfreqs[i]["der"] = freq;
				}
			}
			it++;
		}
		std::cout << i << std::endl;
		std::cout << "anc = " << Cfreqs[i]["anc"] << std::endl;
		std::cout << "der = " << Cfreqs[i]["der"] << std::endl;
	}
}

void fnStats::calcDFreqs(fnFiles &f)
{
	for(unsigned int i=0; i<Dfreqs.size(); i++)
	{
		int total = 0;
		std::unordered_map<std::string,int> l = f.getOutgroupLocus(i);
		std::unordered_map<std::string,int>::iterator it = l.begin();
		while(it != l.end())
		{
			total += it->second;
			it++;
		}

		it = l.begin();
		while(it != l.end())
		{
			double freq = it->second/(double)total;
			if(l.size() == 1){
				if(it->first == ancestral[i])
				{
					Dfreqs[i]["anc"] = freq;
					Dfreqs[i]["der"] = 0.0;
				}
				else if(it->first != ancestral[i])
				{
					Dfreqs[i]["anc"] = 0.0;
					Dfreqs[i]["der"] = freq;
				}
			}
			else if(l.size() == 2)
			{
				if(it->first == ancestral[i])
				{
					Dfreqs[i]["anc"] = freq;
				}
				else
				{
					Dfreqs[i]["der"] = freq;
				}
			}
			it++;
		}
		std::cout << i << std::endl;
		std::cout << "anc = " << Dfreqs[i]["anc"] << std::endl;
		std::cout << "der = " << Dfreqs[i]["der"] << std::endl;
	}
}
