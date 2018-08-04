#include <boost/program_options.hpp> //handling of command line options

#include "fnFiles.h"
#include "fnStats.h"

#include <iostream>
#include <random>
#include <string>

namespace opt = boost::program_options;

//void parseComLine(int argc, char **argv, std::string &infile, std::string &popmap, std::string &outgroup, std::string &abcd, int &vsize);
void parseComLine(int argc, char **argv, std::string &infile, std::string &popmap, std::string &abcd, int &vsize, bool &two, bool &three, bool &four, int &bootstrap);

int main(int argc, char** argv) {

	std::string infile;
	std::string popmap;
	//std::string outgroup;
	std::string abcd;
	int vsize;
	bool two = false;
	bool three = false;
	bool four = false;
	int bootstrap = 0;

	//parseComLine(argc,argv,infile,popmap,outgroup,abcd,vsize);
	parseComLine(argc,argv,infile,popmap,abcd,vsize,two,three,four,bootstrap);

	//std::cout << "Hello World!" << std::endl;

	//start random number generator
	std::default_random_engine generator;
	unsigned int seed = time(0);
	generator.seed(seed);


	//fnFiles f(infile, popmap, outgroup, abcd, vsize);
	fnFiles f(infile, popmap, abcd, vsize);
	f.readfiles();
	if(two == true)
	{
		f.checkF2(); // check to see if all taxa are present in input
	}
	else if (three == true)
	{
		f.checkF3(); // check to see if all taxa are present in input
	}
	else if(four == true)
	{
		f.checkF4(); // check to see if all taxa are present in input
	}
	else
	{
		std::cerr << "No F statistic option was selected." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	fnStats fs(f.getLength(), f.ABCDmap);
	fs.findAncestral(f);
	fs.calcAllFreqs(f, f.ABCDmap);

	if(two == true)
	{
		fs.calcF2(f);
	}
	else if (three == true)
	{
		fs.calcF3(f);
	}
	else if(four == true)
	{
		fs.calcF4(f);
	}
	else
	{
		std::cerr << "No F statistic option was selected." << std::endl;
		exit(EXIT_FAILURE);
	}
	// do bootstrapping
	for(int i=0; i<bootstrap; i++)
	{
		std::vector<int> bootv;
		bootv.resize(f.getLength());
		std::uniform_int_distribution<int> uniform(0,f.getLength()-1);
		for(unsigned int i=0; i<f.getLength(); i++){
			bootv[i] = uniform(generator);
		}
		if(two == true)
		{
			fnStats fsboot(fs,bootv,f.getLength(),f.ABCDmap);
			fnFiles fboot(f,bootv,f.ABCDmap);
			fsboot.calcF2(fboot);
		}
		else if (three == true)
		{
			fnStats fsboot(fs,bootv,f.getLength(),f.ABCDmap);
			fnFiles fboot(f,bootv,f.ABCDmap);
			fsboot.calcF3(fboot);
		}
		else if(four == true)
		{
			fnStats fsboot(fs,bootv,f.getLength(),f.ABCDmap);
			fnFiles fboot(f,bootv,f.ABCDmap);
			fsboot.calcF4(fboot);
		}
		else
		{
			std::cerr << "No F statistic option was selected." << std::endl;
			exit(EXIT_FAILURE);
		}

	}

	return 0;
}

//void parseComLine(int argc, char **argv, std::string &infile, std::string &popmap, std::string &outgroup, std::string &abcd, int &vsize)
void parseComLine(int argc, char **argv, std::string &infile, std::string &popmap, std::string &abcd, int &vsize, bool &two, bool &three, bool &four, int &bootstrap)
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
		("boot,b", opt::value<int>(&bootstrap)->required(), "specify the number of bootstrap replicates to perform.")
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
