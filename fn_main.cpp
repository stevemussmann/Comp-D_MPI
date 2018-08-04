#include <boost/program_options.hpp> //handling of command line options

#include "fnFiles.h"
#include "fnStats.h"

#include <iostream>
#include <string>

namespace opt = boost::program_options;

//void parseComLine(int argc, char **argv, std::string &infile, std::string &popmap, std::string &outgroup, std::string &abcd, int &vsize);
void parseComLine(int argc, char **argv, std::string &infile, std::string &popmap, std::string &abcd, int &vsize, bool &two, bool &three, bool &four);

int main(int argc, char** argv) {

	std::string infile;
	std::string popmap;
	//std::string outgroup;
	std::string abcd;
	int vsize;
	bool two = false;
	bool three = false;
	bool four = false;

	//parseComLine(argc,argv,infile,popmap,outgroup,abcd,vsize);
	parseComLine(argc,argv,infile,popmap,abcd,vsize,two,three,four);

	//std::cout << "Hello World!" << std::endl;

	//fnFiles f(infile, popmap, outgroup, abcd, vsize);
	if(two == true)
	{
		
	}
	else if (three == true)
	{

	}
	else if(four == true)
	{
		fnFiles f(infile, popmap, abcd, vsize);
		f.readfiles();
		f.checkF4(); // check to see if all taxa are present in input

		fnStats fs(f.getLength(), f.ABCDmap);
		fs.findAncestral(f);
		fs.calcAllFreqs(f, f.ABCDmap);
		fs.calcFstats(f);
	}
	else
	{
		std::cerr << "No F statistic option was selected." << std::endl;
		exit(EXIT_FAILURE);
	}

	return 0;
}

//void parseComLine(int argc, char **argv, std::string &infile, std::string &popmap, std::string &outgroup, std::string &abcd, int &vsize)
void parseComLine(int argc, char **argv, std::string &infile, std::string &popmap, std::string &abcd, int &vsize, bool &two, bool &three, bool &four)
{
	opt::options_description desc("--- Option Descriptions ---");
	desc.add_options()
		("help,h", "Prints this help message.")
		("infile,i", opt::value<std::string>(&infile)->required(), "Specifies the input alleles file name.")
		("popmap,p", opt::value<std::string>(&popmap)->required(), "Specifies the population map.")
		("samples,n", opt::value<int>(&vsize)->required(), "Specifies the number of samples in your input file.")
		("abcd,a", opt::value<std::string>(&abcd)->required(), "Specify the file containing a map of your four populations.")
		("two,2", opt::bool_switch(&two), "Use this option to run only the F2 statistic.")
		("three,3", opt::bool_switch(&three), "Use this option to run the F2 and F3 statistics.")
		("four,4", opt::bool_switch(&four), "Use this option to run the F2, F3, and F4 statistics.")
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
	if(two == false && three == false && four == false)
	{
		std::cerr << "You must select one option to calculate F statistics." << std::endl;
		exit(EXIT_FAILURE);
	}
	if(two == true && three == true && four == true)
	{
		std::cerr << "You must select only one option to calculate F statistics." << std::endl;
		exit(EXIT_FAILURE);
	}

	if((two == true && three == true) || (three == true && four == true) || (two == true && four == true))
	{
		std::cerr << "You must select only one option to calculate F statistics." << std::endl;
		exit(EXIT_FAILURE);
	}
	
}
