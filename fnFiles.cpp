#include "fnFiles.h"

//#include <algorithm>
//#include <cstdlib>
#include <fstream>
#include <functional> //for reverse order of keys in map
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

//fnFiles::fnFiles(std::string i, std::string p, std::string o, std::string abcd, int vectorsize)
fnFiles::fnFiles(std::string i, std::string p, std::string abcd, int vectorsize)
{
	infile = i;
	popfile = p;
	ABCDfile = abcd;
	readABCDfile(); //must read ABCDfile in the constructor to get taxa invovled
	readPopfile();
	for(std::unordered_map<std::string,std::string>::iterator it = ABCDmap.begin(); it != ABCDmap.end(); it++)
	{
		data[it->second].resize(vectorsize);
	}

}

std::unordered_map <std::string,int> fnFiles::getOutgroupLocus(int i)
{
	return data["D"][i];
}

std::unordered_map <std::string,int> fnFiles::getALocus(int i)
{
	return data["A"][i];
}

std::unordered_map <std::string,int> fnFiles::getBLocus(int i)
{
	return data["B"][i];
}

std::unordered_map <std::string,int> fnFiles::getCLocus(int i)
{
	return data["C"][i];
}

int fnFiles::getLength()
{
	if(data["A"].size() == data["B"].size() && data["B"].size() == data["C"].size() && data["C"].size() == data["D"].size())
	{
		return data["A"].size();
	}
	else
	{
		std::cerr << "Vectors holding data for taxa A,B,C, and D are different lengths." << std::endl;
		exit(EXIT_FAILURE);
	}
}

void fnFiles::readfiles()
{
	std::cout << "Reading files" << std::endl;
	readPhylip();
	std::cout << "blacklisting loci" << std::endl;
	blacklist();
}

void fnFiles::readPhylip()
{
	std::ifstream myfile(infile.c_str()); //convert file to stream
	int locnumber;
	int counter=0;

	if(myfile.is_open())
	{
		std::string line;
		while(getline(myfile, line))
		{
			std::vector<std::string> tokens; //vector to hold split line
			std::istringstream iss(line); //convert to stream
			copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), back_inserter(tokens)); //put stream into vector
			
			if(counter ==0)
			{
				std::stringstream num(tokens[1]);
				num >> locnumber;
			}
			else
			{
				std::unordered_map<std::string,std::string>::const_iterator got = popmap.find(tokens[0]);

				if(got == popmap.end() )
				{
					std::cout << "WARNING: sample " << tokens[0] << " not found in popmap and will not be used in calculations" << std::endl;
					std::cout << "Verify that this is OK before interpreting your results." << std::endl << std::endl;
				}
				else
				{
					std::unordered_map<std::string,std::string>::const_iterator got2 = ABCDmap.find(popmap[tokens[0]]);

					if(got2 == ABCDmap.end() )
					{
						std::cout << "Sample " << tokens[0] << " from population " << popmap[tokens[0]] << " is being ignored." << std::endl << std::endl;
					}
					else
					{

						if(tokens[1].size() != (unsigned int)locnumber){
							std::cerr << "Length of Phylip sequence is " << tokens[1].size() << " but input length was " << locnumber << std::endl;
							exit(EXIT_FAILURE);
						}
						for(int i=0; i<locnumber; i++) //put species name into locusfile object
						{
							//std::cout << "Putting data into data structure" << std::endl;
	
							std::stringstream ss;
							std::string tempstring;
							ss << tokens[1][i];
							ss >> tempstring;
	
							if(tempstring == "M" || tempstring == "R" || tempstring == "W" || tempstring == "S" || tempstring == "Y" || tempstring == "K")
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

								// new data structure
								data[ABCDmap[popmap[tokens[0]]]][i][allele0]+=1;
								data[ABCDmap[popmap[tokens[0]]]][i][allele1]+=1;
		
							}
							else if( tempstring == "A" || tempstring == "C" || tempstring == "G" || tempstring == "T")
							{
								// new data structure
								data[ABCDmap[popmap[tokens[0]]]][i][tempstring]+=2;

							}
						}
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
		std::cout.flush();
		exit(EXIT_FAILURE);
	}
}

void fnFiles::readPopfile()
{
	std::ifstream myfile(popfile.c_str());
	if(myfile.is_open())
	{
		std::string line;
		while(getline(myfile,line))
		{
			//std::cout << line << std::endl;
			std::vector<std::string> tokens;
			std::istringstream iss(line);
			copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), back_inserter(tokens));
			
			popmap[tokens[0]] = tokens[1];
		}
	}
	else
	{
		std::cerr << "ERROR: File " << popfile << " not found." << std::endl;
		exit(EXIT_FAILURE);
	}
	myfile.close();
}

void fnFiles::readABCDfile()
{
	std::ifstream myfile(ABCDfile.c_str());
	if(myfile.is_open())
	{
		std::string line;
		while(getline(myfile,line))
		{
			//std::cout << line << std::endl;
			std::vector<std::string> tokens;
			std::istringstream iss(line);
			copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), back_inserter(tokens));
			
			ABCDmap[tokens[0]] = tokens[1];
		}
	}
	else
	{
		std::cerr << "ERROR: File " << ABCDfile << " not found." << std::endl;
		exit(EXIT_FAILURE);
	}
	myfile.close();
}

std::string fnFiles::iupac(std::string ambig)
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

void fnFiles::blacklist()
{
	std::map<int,int, std::greater<int> > bl; //list of blacklisted loci

	for(unsigned int i=0; i<data["A"].size(); i++)
	{
		for(std::unordered_map<std::string,std::string>::iterator it = ABCDmap.begin(); it != ABCDmap.end(); it++)
		{
			if(data[it->second][i].size() < 1 || data[it->second][i].size() > 2)
			{
				bl[i]++;
			}
		}
	}

	//remove loci from vectors
	std::map<int,int>::iterator it = bl.begin();
	while(it != bl.end())
	{
		for(std::unordered_map<std::string,std::string>::iterator itt = ABCDmap.begin(); itt != ABCDmap.end(); itt++)
		{
			data[itt->second].erase(data[itt->second].begin()+it->first);
		}
		it++;
	}
}
