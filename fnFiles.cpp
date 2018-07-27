#include "fnFiles.h"

//#include <algorithm>
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

fnFiles::fnFiles(std::string i, std::string p, std::string o, std::string abcd, int vectorsize){
	infile = i;
	popfile = p;
	ABCDfile = abcd;
	outgroup = o;
	A.resize(vectorsize);
	B.resize(vectorsize);
	C.resize(vectorsize);
	D.resize(vectorsize);

}

void fnFiles::readfiles()
{
	std::cout << "Reading files" << std::endl;
	readABCDfile();
	readPopfile();
	readPhylip();
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
							//species[i].push_back(tokens[0]);
	
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
		
								if(ABCDmap[popmap[tokens[0]]] == "A")
								{
									A[i][allele0]+=1;
									A[i][allele1]+=1;
								}
								else if(ABCDmap[popmap[tokens[0]]] == "B")
								{
									B[i][allele0]+=1;
									B[i][allele1]+=1;
								}
								else if(ABCDmap[popmap[tokens[0]]] == "C")
								{
									C[i][allele0]+=1;
									C[i][allele1]+=1;
								}
								else if(ABCDmap[popmap[tokens[0]]] == "D")
								{
									D[i][allele0]+=1;
									D[i][allele1]+=1;
								}
								else
								{
									std::cout << "This code should be unreachable" << std::endl;
									exit(EXIT_FAILURE);
								}
							}
							else if( tempstring == "A" || tempstring == "C" || tempstring == "G" || tempstring == "T")
							{
								if(ABCDmap[popmap[tokens[0]]] == "A")
								{
									A[i][tempstring]+=2;
								}
								else if(ABCDmap[popmap[tokens[0]]] == "B")
								{
									B[i][tempstring]+=2;
								}
								else if(ABCDmap[popmap[tokens[0]]] == "C")
								{
									C[i][tempstring]+=2;
								}
								else if(ABCDmap[popmap[tokens[0]]] == "D")
								{
									D[i][tempstring]+=2;
								}
								else
								{
									std::cout << "This code should be unreachable" << std::endl;
									exit(EXIT_FAILURE);
								}
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
