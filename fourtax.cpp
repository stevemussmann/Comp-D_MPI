/* 
 * File:   locus.cpp
 * Author: Steve
 * 
 * Created on July 7, 2015, 9:03 PM
 */

#include "fourtax.h"
#include "locusfile.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <random>

fourtax::fourtax(int locusnum, int taxnum) {
	number.resize(locusnum);
	pattern.resize(locusnum);
	names.resize(locusnum);
	seq.resize(locusnum);
	major.resize(locusnum);
	freq.resize(locusnum); //add a vector of vector<double> to hold allele frequencies for each individual
	for(int i=0; i<locusnum; i++)
	{
		names[i].resize(taxnum);
		seq[i].resize(taxnum);
		major[i].resize(taxnum);
		freq[i].resize(taxnum);
	}
}

fourtax::fourtax(const fourtax& orig) {
}

std::string fourtax::getPattern(int i)
{
    return pattern.at(i);
}

double fourtax::getFreq(int locus, int taxon)
{
	return freq[locus].at(taxon);
}

void fourtax::calculatePattern(int locus, int ntaxa, locusfile &file, int my_rank)
{
	for(int i=ntaxa-1; i>0; i--)
	{
		if(seq[locus].at(i) == major[locus].at(0)) //compare to major allele for outgroup
		{
			pattern[locus].append("A");
		}
		else
		{
			pattern[locus].append("B");
		}
	}
	pattern[locus].append("A");
}

void fourtax::populateDtest(std::vector<int> &keep, locusfile &file, std::unordered_map <std::string,int> &indlist, 
			    std::default_random_engine &generator, int my_rank, int ntaxa)
{
	std::binomial_distribution<int> binomial(1,0.5);
    
	std::unordered_map<std::string,int>::iterator findit;
	for(unsigned int i=0; i<keep.size(); i++) //for each kept locus
	{
		for(int j=0; j<file.GetSize(keep[i]); j++) //for each individual
		{
			std::string name = file.GetName(keep[i],j);
			std::string gene;
			if(file.GetSeq(keep[i],j*2) != file.GetSeq(keep[i],(j*2)+1))
			{
				int number = binomial(generator);
				if(number == 0)
				{
					gene = file.GetSeq(keep[i],j*2);
				}
				else
				{
					gene = file.GetSeq(keep[i],(j*2)+1);
				}
			}
			else
			{
				gene = file.GetSeq(keep[i],j*2);
			}
			//cout << gene << endl;
			findit = indlist.find(name);
			if(findit != indlist.end())
			{
				//cout << file[keep[i]].GetName(j) << " is " << findit->second << endl;
				names[i].at(findit->second) = name;
				seq[i].at(findit->second) = gene;
				//dtest[i].setSeq(gene, findit->second);               
			}
		}
		
		//add taxon major alleles
		for(int j=0; j<ntaxa; j++)
		{
			major[i].at(j) = file.getMajorAllele(keep[i], j, my_rank);
		}
		
	}
}

void fourtax::populateDtest(std::vector<int> &keep, locusfile &file, std::unordered_map<std::string,int> &indlist, int my_rank, int ntaxa)
{
	std::unordered_map<std::string,int>::iterator findit;
	for(unsigned int i=0; i<keep.size(); i++)//for each kept locus
	{
		for(int j=0; j<ntaxa; j++)
		{
			major[i].at(j) = file.getMajorAllele(keep[i], j, my_rank); //find major allele for each taxon at each locus
		}
		
		for(int j=0; j<file.GetSize(keep[i]); j++)//for each individual
		{
			std::string name = file.GetName(keep[i],j); //get sample name
			//find major allele from population information
			
			findit = indlist.find(name);
			if(findit != indlist.end())
			{
				names[i].at(findit->second) = name;
				//do something to determine if heterozygote, and if one allele or other matches outgroup
				std::string allele1 = file.GetSeq(keep[i],j*2);
				std::string allele2 = file.GetSeq(keep[i],(j*2)+1);
				
				if(allele1 == major[i].at(0) && allele2 == major[i].at(0))
				{
					freq[i].at(findit->second) = 0.0; //is a homozygote match for ancestral allele
				}
				else if(allele1 != major[i].at(0) && allele2 != major[i].at(0))
				{
					freq[i].at(findit->second) = 1.0; //is a homozygote match for derived allele
				}
				else
				{
					freq[i].at(findit->second) = 0.5; //is heterozygote
				}
			}
		}
	}
}

fourtax::~fourtax() {
}
