/* 
 * File:   locusfile.cpp
 * Author: Steve
 * 
 * Created on July 7, 2015, 9:03 PM
 */

#include "locusfile.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <cassert>

#include "mpi.h"

locusfile::locusfile(int vectorsize)
{
	species.resize(vectorsize);
	alleles.resize(vectorsize);
	freqs.resize(vectorsize);
}

void locusfile::AddData(int i, std::string name, std::string sequence, int allele)
{
	if(allele == 0)
	{
		species[i].push_back(name);
	}
	//put the sequence onto the vector
	alleles[i].push_back(sequence);
}

void locusfile::AddData(int i, std::string name, std::string allele1, std::string allele2, std::vector<std::unordered_map<std::string,int> > freq)
{
	species[i].push_back(name);
	alleles[i].push_back(allele1);
	alleles[i].push_back(allele2);
	freqs[i]=freq;
}

int locusfile::GetSize(int i)
{
	int length = species[i].size();
	return length;
}

int locusfile::GetSeqSize(int i)
{
	int length = alleles[i].size();
	return length;
}

std::string locusfile::GetSeq(int i, int j)
{
	std::string seq = alleles[i].at(j);
	return seq;
}

std::string locusfile::GetName(int i, int j)
{
	std::string name = species[i].at(j);
	return name;
}

//for reading pyrad .alleles files
void locusfile::readInput(std::string infile, int locnumber)
{
	//read the input file
	std::ifstream myfile(infile.c_str());
	if( myfile.is_open()) //check if file is open
	{
		int loc=0; //initialize locus counter
		std::string sep="|"; // initialize locus separator
		std::string line;
		while(getline(myfile,line))
		{
			//std::string line; //declare temporary string for line
			//getline(myfile,line); //capture each line
			size_t found=line.find(sep); //check for locus separator
			
			if(found!=std::string::npos) //if locus separator present
			{
				loc++; //increment locus counter
				if( loc > locnumber) //error message to print if the user defined number of alleles is less than the actual number
				{
					std::cerr << "Input pyrad alleles file contains more than " << locnumber << " alleles" << std::endl;
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				std::string a, name, seq, allele;//a is temporary

				line.erase(line.begin()); //remove first character

				std::stringstream ss(line);
				ss >> a >> seq; //name and allele go to a, sequence goes to seq

				//strip characters
				a.replace(a.end()-2, a.end()-1, " "); //replace underscore with space
				
				std::stringstream ss2(a); //new stringstream
				
				ss2 >> name >> allele; //split to name and allele
				
				int iallele; //stores allele number as int
				std::istringstream(allele) >> iallele; //convert allele to int
				
				if(iallele == 0)
				{
					species[loc].push_back(name);
				}
				//put the sequence onto the vector
				alleles[loc].push_back(seq);
			}
		}
		myfile.close(); //close file
	}
	else //kill program if input file cannot be read
	{
		std::cerr << "Unable to open " << infile << std::endl;
		MPI_Finalize();
		std::cout.flush();
		exit(EXIT_FAILURE);
	}
}

//for reading structure files
void locusfile::readInput(std::string infile, int locnumber, int offset)
{
	std::ifstream myfile(infile.c_str()); //convert file to stream
	offset=offset+1;
	int counter=0; //counter for line number of file
	
	if( myfile.is_open())
	{
		std::string line;//string to temporarily hold line
		while(getline(myfile, line))
		{
			std::vector<std::string> tokens; //vector to temporarily hold split line
			std::istringstream iss(line); //convert to stream
			copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), back_inserter(tokens)); //put stream into vector
			
			if(counter%2 == 0) //on even numbered lines, put the species name onto the species vector
			{
				for(int i=0; i<locnumber; i++)
				{
					species[i].push_back(tokens[0]);
				}
			}
			for(int i=0; i<locnumber; i++)
			{
				alleles[i].push_back(tokens[i+offset]);
			}
			counter++;
		}
		myfile.close(); //close input file
	}
	else
	{
		std::cerr << "Unable to open " << infile << std::endl;
		MPI_Finalize();
		std::cout.flush();
		exit(EXIT_FAILURE);
	}
}

void locusfile::readInput(std::string infile, int locnumber, bool phylip)
{
	std::ifstream myfile(infile.c_str()); //convert file to stream
	int counter=0;
	
	if(myfile.is_open())
	{
		std::string line;
		while(getline(myfile, line))
		{
			if(counter != 0)
			{
				std::vector<std::string> tokens; //vector to hold split line
				std::istringstream iss(line); //convert to stream
				copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), back_inserter(tokens)); //put stream into vector
				//std::cout << tokens[1].size() << std::endl; //print size of tokens vector to screen
				if(tokens[1].size() != (unsigned int)locnumber){
					std::cerr << "Length of Phylip sequence is " << tokens[1].size() << " but input length was " << locnumber << std::endl;
					exit(EXIT_FAILURE);
				}
				for(int i=0; i<locnumber; i++) //put species name into locusfile object
				{
					species[i].push_back(tokens[0]);
					
					std::stringstream ss;
					std::string tempstring;
					ss << tokens[1][i];
					ss >> tempstring;
					
					if(tempstring == "M" || tempstring == "R" || tempstring == "W" || 
						tempstring == "S" || tempstring == "Y" || tempstring == "K")
					{
						std::string bases = iupac(tempstring);
						
						std::stringstream as0;
						std::stringstream as1;
						std::string allele0;
						std::string allele1;
						as0 << bases[0];
						as1 << bases[1];
						as0 >> allele0;
						as1 >> allele1;
						
						alleles[i].push_back(allele0);
						alleles[i].push_back(allele1);
					}
					else
					{
						alleles[i].push_back(tempstring);
						alleles[i].push_back(tempstring);
					}
				}
			}
			counter++;
		}
		myfile.close();
	}
	else
	{
		std::cerr << "Unable to open " << infile << std::endl;
		MPI_Finalize();
		std::cout.flush();
		exit(EXIT_FAILURE);
	}
}

void locusfile::removeN(int locnumber, bool Nremoveflag, bool gapignoreflag)
{
	char findN='N';
	char findDash='-';
	
	//find all occurrences of N at a locus
	for(int i=0; i<locnumber; i++) //iterate over each locus
	{
		std::vector<int> allelevector; //vector to hold positions to be removed from each allele
		for(unsigned int j=0; j<alleles[i].size(); j++) //iterate over each allele at each locus
		{
			//run findLocation subroutine on each allele
			std::vector<int> tempvector = findLocation(alleles[i].at(j), findN, findDash, Nremoveflag, gapignoreflag);
			//append tempvector to allelevector
			allelevector.insert( allelevector.end(), tempvector.begin(), tempvector.end() );
		}
		// sort the vector of character positions 
		std::sort(allelevector.begin(), allelevector.end());
		
		//find just the unique character position numbers
		std::vector<int>::iterator it;
		it=unique(allelevector.begin(), allelevector.end());
		
		//resize the allelevector to have only the unique numbers
		allelevector.resize(distance(allelevector.begin(), it));

		//replace the unique positions at the locus with nothing
		for(unsigned int j=0; j<alleles[i].size(); j++) //iterate over each allele at each locus
		{
			std::string tempstring=alleles[i].at(j);
			
			//for each allele at the locus, remove the positions
			if(allelevector.size() != 0)
			{
				//count backwards through the array to avoid shifting positions
				for(int pos=allelevector.size()-1; pos>=0; pos--)
				{
					//erase the positions containing unwanted characters
					tempstring.erase(allelevector[pos], 1);
				}
			}
			alleles[i].at(j) = tempstring; // replace the original allele with the corrected allele
		}
	}
}

std::vector<int> locusfile::findLocation(std::string sequence, char findN, char findDash, bool Nremoveflag, bool gapignoreflag)
{
	std::vector<int> NLocations;
	for(unsigned int i=0; i<sequence.size(); i++)
	{
		if(Nremoveflag==true)
		{
			if(sequence[i] == findN)
			{
				NLocations.push_back(i);
			}
		}
		if(gapignoreflag==true)
		{
			if(sequence[i] == findDash)
			{
				NLocations.push_back(i);
			}
		}
	    
	}
	return NLocations;
}

std::vector<int> locusfile::findInformative(locusfile &current, int locnumber, std::unordered_map <std::string,int> indlist, int ntaxa, bool hetIgnore)
{
	std::vector<int> keep;

	//create a vector of locus objects with data from the tiptaxa
	for(int i=0; i<locnumber; i++)
	{
		for(unsigned int j=0; j<species[i].size(); j++)
		{
			std::unordered_map<std::string, int>::const_iterator got=indlist.find(species[i].at(j));
			if(got!=indlist.end())
			{   
				if(alleles[i].at(j*2) != "-9" && alleles[i].at((j*2)+1) != "-9" && 
					alleles[i].at(j*2) != "N" && alleles[i].at((j*2)+1) != "N" && 
					alleles[i].at(j*2) != "-" && alleles[i].at((j*2)+1) != "-")
				{
					current.AddData(i, species[i].at(j), alleles[i].at(j*2), alleles[i].at((j*2)+1), freqs[i]);
				}
			}
		}

		if(current.GetSize(i) == ntaxa)
		{
			std::unordered_map <std::string,int> seq_count; //map for finding unique allele copies
			bool heterozygote = false;
			for(int j=0; j<current.GetSeqSize(i); j+=2)
			{
				if(current.GetSeq(i,j) == current.GetSeq(i,j+1)) //if individual is homozygote
				{
					seq_count[current.GetSeq(i,j)]++; //add to the map
				}
				else //if heterozygote
				{
					heterozygote = true;
					seq_count[current.GetSeq(i,j)]++; //add allele 1 to the map
					seq_count[current.GetSeq(i,j+1)]++; //add allele 2 to the map                    
				}
			}
            
			if(hetIgnore==true)
			{
				if(seq_count.size() == 2 && heterozygote == false) //if locus is biallelic and all inds are homozygotes
				{
					keep.push_back(i); //add locus number to keep vector
				}
			}
			else
			{
				if(seq_count.size() == 2)
				{
					keep.push_back(i);
				}
			}
		}
	}
	return keep;
}

std::string locusfile::iupac(std::string ambig)
{
	std::unordered_map <std::string, std::string> map;
	map["M"] = "AC";
	map["R"] = "AG";
	map["W"] = "AT";
	map["S"] = "CG";
	map["Y"] = "CT";
	map["K"] = "GT";
	
	return map[ambig];
}

void locusfile::calcFreq(int locnumber, std::vector<std::vector<std::string>> &taxa)
{
	for(unsigned int i=0; i<taxa.size(); i++){ //for each group of taxa (o, p1, p2, etc...)
		std::unordered_map<std::string,int> taxList; //make a list of taxa
		for(unsigned int j=0; j<taxa[i].size(); j++)
		{
			taxList[taxa[i][j]] = 1;
		}
		
		for(int k = 0; k<locnumber; k++) //for each locus
		{
			std::unordered_map<std::string,int>alleleCount; //keep track of alleles at a locus for a group of taxa
			for(unsigned int l=0; l<species[k].size(); l++) //for each individual, check to see if it is in the list
			{
				std::unordered_map<std::string,int>::const_iterator got=taxList.find(species[k].at(l));
				if(got!=taxList.end())
				{
					if(alleles[k].at(l*2) != "-9" && alleles[k].at(l*2) != "N" && alleles[k].at(l*2) != "-")
					{
						alleleCount[alleles[k].at((l*2))]++;
						alleleCount[alleles[k].at((l*2)+1)]++;
					}
				}
			}
			freqs[k].push_back(alleleCount);
		}
	}
}

std::string locusfile::getMajorAllele(int locus, int taxon, int my_rank)
{
	std::string major="";
	
	int max = 0;
	
	if(freqs[locus][taxon].size() == 1){
		for(auto it = freqs[locus][taxon].begin(); it != freqs[locus][taxon].end(); it++)
		{
			major = it->first;
		}
	}
	else
	{
		std::unordered_map<int,int> counts;
		for(auto it = freqs[locus][taxon].begin(); it != freqs[locus][taxon].end(); it++)
		{
			counts[it->second]++;
		}
		
		if(counts.size() == 1)
		{
			for(auto it = freqs[locus][taxon].begin(); it != freqs[locus][taxon].end(); it++)
			{
				major = it->first;
			}
		}
		else
		{
			for(auto it = freqs[locus][taxon].begin(); it != freqs[locus][taxon].end(); it++)
			{
				if(it->second > max)
				{
					major = it->first;
					max = it->second;
				}
			}
		}
	}
	assert(major != "");
	
	return major;
}

locusfile::~locusfile() {
}