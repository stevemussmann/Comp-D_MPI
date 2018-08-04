#include "fnFiles.h"
#include "fnStats.h"

//#include <algorithm>
//#include <cstdlib>
//#include <fstream>
#include <iomanip>
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

/*
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
*/

void fnStats::findAncestral(fnFiles &f)
{
	for(unsigned int i=0; i<ancestral.size(); i++)
	{
		std::unordered_map<std::string,int> l = f.getDLocus(i); //get locus for outgroup
		std::unordered_map<std::string,int>::iterator it = l.begin(); //start iterator
		ancestral[i] = it->first; //if only one allele, set it as ancestral state.
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
		//std::cout << i << std::endl;
		//std::cout << "anc = " << Afreqs[i]["anc"] << std::endl;
		//std::cout << "der = " << Afreqs[i]["der"] << std::endl;
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
		//std::cout << i << std::endl;
		//std::cout << "anc = " << Bfreqs[i]["anc"] << std::endl;
		//std::cout << "der = " << Bfreqs[i]["der"] << std::endl;
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
		//std::cout << i << std::endl;
		//std::cout << "anc = " << Cfreqs[i]["anc"] << std::endl;
		//std::cout << "der = " << Cfreqs[i]["der"] << std::endl;
	}
}

void fnStats::calcDFreqs(fnFiles &f)
{
	for(unsigned int i=0; i<Dfreqs.size(); i++)
	{
		int total = 0;
		std::unordered_map<std::string,int> l = f.getDLocus(i);
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
		//std::cout << i << std::endl;
		//std::cout << "anc = " << Dfreqs[i]["anc"] << std::endl;
		//std::cout << "der = " << Dfreqs[i]["der"] << std::endl;
	}
}

void fnStats::calcFstats(fnFiles &f)
{
	double f2total = 0.0;
	double f3total = 0.0;
	double f4total = 0.0;

	for(unsigned int i=0; i<ancestral.size(); i++)
	{
		//get locus counts	
		std::unordered_map<std::string,int> al = f.getALocus(i);
		std::unordered_map<std::string,int> bl = f.getBLocus(i);
		std::unordered_map<std::string,int> cl = f.getCLocus(i);
		std::unordered_map<std::string,int> dl = f.getDLocus(i);

		//initialize iterators
		std::unordered_map<std::string,int>::iterator ait = al.begin();
		std::unordered_map<std::string,int>::iterator bit = bl.begin();
		std::unordered_map<std::string,int>::iterator cit = cl.begin();
		std::unordered_map<std::string,int>::iterator dit = dl.begin();

		std::vector<int> avec;
		std::vector<int> bvec;
		std::vector<int> cvec;
		std::vector<int> dvec;

		int atotal=0;
		while(ait != al.end())
		{
			atotal += ait->second;
			avec.push_back(ait->second);
			ait++;
		}

		int btotal=0;
		while(bit != bl.end())
		{
			btotal += bit->second;
			bvec.push_back(bit->second);
			bit++;
		}

		int ctotal=0;
		while(cit != cl.end())
		{
			ctotal += cit->second;
			cvec.push_back(cit->second);
			cit++;
		}

		int dtotal=0;
		while(dit != dl.end())
		{
			dtotal += dit->second;
			dvec.push_back(dit->second);
			dit++;
		}

		double ahz = hz(avec, atotal);
		double bhz = hz(bvec, btotal);
		double chz = hz(cvec, ctotal);
		double dhz = hz(dvec, dtotal);

		double f2i = f2(Afreqs[i]["anc"], Bfreqs[i]["anc"], ahz, bhz, atotal, btotal);
		f2total+=f2i;
		double f3i = f3(Afreqs[i]["anc"], Bfreqs[i]["anc"], Cfreqs[i]["anc"], chz, ctotal);
		f3total+=f3i;
		double f4i = f4(Afreqs[i]["anc"], Bfreqs[i]["anc"], Cfreqs[i]["anc"], Dfreqs[i]["anc"]);
		f4total+=f4i;
		double fsti = fst(Afreqs[i]["der"], Bfreqs[i]["der"]);

		//std::cout << std::fixed;
		std::cout << "len(a) = " << avec.size() << std::endl;
		std::cout << "len(b) = " << bvec.size() << std::endl;
		std::cout << "len(c) = " << cvec.size() << std::endl;
		std::cout << "ahz = " << ahz << std::endl;
		std::cout << "bhz = " << bhz << std::endl;
		std::cout << "chz = " << chz << std::endl;
		std::cout << "f2 = " << f2i << std::endl;
		std::cout << "f3 = " << f3i << std::endl;
		//std::cout << std::setprecision(10) << "f4 = " << f4i << std::endl;
		std::cout << "f4 = " << f4i << std::endl;
		std::cout << "fst = " << fsti << std::endl;
		std::cout << std::endl;
		
	}
	std::cout << std::endl << std::endl;
	double f2avg = f2total/(double)ancestral.size();
	double f3avg = f3total/(double)ancestral.size();
	double f4avg = f4total/(double)ancestral.size();
	std::cout << "f2avg = " << f2avg << std::endl;
	std::cout << "f3avg = " << f3avg << std::endl;
	std::cout << "f4avg = " << f4avg << std::endl;
}

double fnStats::hz(std::vector<int> &v, int t)
{
	if(v.size() == 1){
		return 0.0;
	}
	else
	{
		return ((double)v[0]*(double)v[1])/((double)t*((double)t-1.0));
	}
}


double fnStats::f2(double freqa, double freqb, double ha, double hb, int tota, int totb)
{
	double x = freqa-freqb;
	return (x*x)-(ha/(double)tota)-(hb/(double)totb);
}

double fnStats::f3(double freqa, double freqb, double freqc, double hc, int totc)
{
	double x = freqc-freqa;
	double y = freqc-freqb;
	double xy = (x*y);
	return xy-(hc/(double)totc);
}

double fnStats::f4(double freqa, double freqb, double freqc, double freqd)
{
	double x = freqa-freqb;
	double y = freqc-freqd;
	return x*y;
}

double fnStats::fst(double freqa, double freqb)
{
	double n = (freqa-freqb);
	n = n*n;
	double b = 1.0-freqb;
	double a = 1.0-freqa;
	b = freqa*b;
	a = freqb*a;
	if(a+b == 0)
	{
		return 0.0;
	}
	else
	{
		return n/(a+b);
	}
}
