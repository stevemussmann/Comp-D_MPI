#include <boost/program_options.hpp> //handling of command line options

#include "fnFiles.h"
#include "fnStats.h"

#include <iostream>
#include <string>

namespace opt = boost::program_options;

//void parseComLine(int argc, char **argv, std::string &infile, std::string &popmap, std::string &outgroup, std::string &abcd, int &vsize);
void parseComLine(int argc, char **argv, std::string &infile, std::string &popmap, std::string &abcd, int &vsize);

int main(int argc, char** argv) {

	std::string infile;
	std::string popmap;
	//std::string outgroup;
	std::string abcd;
	int vsize;

	//parseComLine(argc,argv,infile,popmap,outgroup,abcd,vsize);
	parseComLine(argc,argv,infile,popmap,abcd,vsize);

	//std::cout << "Hello World!" << std::endl;

	//fnFiles f(infile, popmap, outgroup, abcd, vsize);
	fnFiles f(infile, popmap, abcd, vsize);
	f.readfiles();

	fnStats fs(f.getLength(), f.ABCDmap);
	fs.findAncestral(f);
	fs.calcAllFreqs(f, f.ABCDmap);
	fs.calcFstats(f);

	return 0;
}

//void parseComLine(int argc, char **argv, std::string &infile, std::string &popmap, std::string &outgroup, std::string &abcd, int &vsize)
void parseComLine(int argc, char **argv, std::string &infile, std::string &popmap, std::string &abcd, int &vsize)
{
	opt::options_description desc("--- Option Descriptions ---");
	desc.add_options()
		("help,h", "Prints this help message.")
		("infile,i", opt::value<std::string>(&infile)->required(), "Specifies the input alleles file name.")
		("popmap,p", opt::value<std::string>(&popmap)->required(), "Specifies the population map.")
		//("outgroup,o", opt::value<std::string>(&outgroup)->required(), "Specifies the outgroup taxon name.")
		("samples,n", opt::value<int>(&vsize)->required(), "Specifies the number of samples in your input file.")
		("abcd,a", opt::value<std::string>(&abcd)->required(), "Specify the file containing a map of your four populations.")
	;

	opt::variables_map vm;
	try
	{
		opt::store(opt::parse_command_line(argc,argv,desc), vm);

		if(vm.count("help"))
		{
			std::cout << "This program was created for calculating f statistics." << std::endl;
			std::cout << desc << std::endl;
			exit(EXIT_FAILURE);
		}
		opt::notify(vm); // throws an error if there are any problems

	}
	catch(opt::required_option& e) //catch errors resulting from required options
	{
		std::cerr << std::endl << "ERROR: " << e.what() << std::endl << std::endl;
		std::cout << desc << std::endl;
		exit(EXIT_FAILURE);
	}
	catch(opt::error& e) //catch other command line errors
	{
		std::cerr << std::endl << "ERROR: " << e.what() << std::endl << std::endl;
		std::cout << desc << std::endl;
		exit(EXIT_FAILURE);
	}
	
}
