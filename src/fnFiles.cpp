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

// copy constructor for bootstrap replication
fnFiles::fnFiles(fnFiles f, std::vector<int> &v, std::unordered_map<std::string,std::string> m)
{
	for(std::unordered_map<std::string,std::string>::iterator it = m.begin(); it != m.end(); it++)
	{
		data[it->second].resize(v.size());
	}
	
	for(unsigned int i=0; i<v.size(); i++)
	{
		for(std::unordered_map<std::string,std::string>::iterator it = m.begin(); it!=m.end(); it++)
		{
			std::unordered_map<std::string,int> la = f.data[it->second][v[i]];
			data[it->second][i] = la;
		}
	}

}

void fnFiles::checkF2()
{
	std::vector<std::string> v = {"A","B","O"};
	checkPops(v);
}

void fnFiles::checkF3()
{
	std::vector<std::string> v = {"A","B","C","O"};
	checkPops(v);
}

void fnFiles::checkF4()
{
	std::vector<std::string> v = {"A","B","C","D","O"};
	checkPops(v);
}

void fnFiles::checkPops(std::vector<std::string> &v)
{
	std::unordered_map<std::string,std::string> temp;
	for(std::unordered_map<std::string,std::string>::iterator it = ABCDmap.begin(); it != ABCDmap.end(); it++)
	{
		temp[it->second] = it->first;
	}

	for(unsigned int i=0; i<v.size(); i++){

		std::unordered_map<std::string,std::string>::const_iterator got = temp.find(v.at(i));
		if(got == temp.end())
		{
			std::cout << "Taxon " << v.at(i) << " was not found in your population file." << std::endl;
			std::cout << "Please check your population file." << std::endl << std::endl;
			exit(EXIT_FAILURE);
		}

	}
}

std::unordered_map <std::string,int> fnFiles::getLocus(std::string s, int i)
{
	return data[s][i];
}

unsigned int fnFiles::getLength()
{
	int counter=0;
	unsigned int now = 0;
	for(std::unordered_map<std::string,std::string>::iterator it = ABCDmap.begin(); it != ABCDmap.end(); it++)
	{
		if(counter == 0)
		{
			now = data[it->second].size();
			//std::cout << now << std::endl;
		}
		else
		{
			if( now != data[it->second].size() )
			{
				std::cerr << "Vectors holding data for taxa A,B,C, and D are different lengths." << std::endl;
				exit(EXIT_FAILURE);
			}
			else
			{
				now = data[it->second].size();
			}
		}
		counter++;
	}

	return now;
}

void fnFiles::readfiles(int vectorsize, bool p, bool s, std::string missing, int offset)
{
	std::cout << "Reading files" << std::endl;
	//std::string missing = "-9";
	if(s==true)
	{
		readStructure(offset,vectorsize,missing);
	}
	else if(p==true)
	{
		readPhylip();
	}
	else
	{
		std::cerr << "Booleans for both structure and phylip file formats are true.  This code block should be unreachable - how did you get here?" << std::endl;
		exit(EXIT_FAILURE);
	}
	std::cout << "blacklisting loci" << std::endl;
	blacklist();
}

void fnFiles::readStructure(int offset, int l, std::string m)
{
	std::ifstream myfile(infile.c_str()); //convert file to stream
	offset=offset+1;
	int counter=0; //counter for line number of file

	if(myfile.is_open())
	{
		std::string line; //string to temporarily hold line
		while(getline(myfile,line))
		{
			std::vector<std::string> tokens; //vector to temporarily hold split line
			std::istringstream iss(line); //convert to stream
			copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), back_inserter(tokens)); //put stream into vector
				
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

					//if(counter%2 == 0) //on even numbered lines, put species name onto the species vector
					for(int i=0; i<l; i++)
					{
						//std::cout << i << std::endl;
						if(tokens[i+offset] != m) //exclude missing data
						{
							data[ABCDmap[popmap[tokens[0]]]][i][tokens[i+offset]]+=1;
						}
					}	
					counter++;
				}
			}
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
