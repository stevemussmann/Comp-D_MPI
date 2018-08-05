#include <boost/program_options.hpp> //handling of command line options

#include "fnFiles.h"
#include "fnStats.h"

#include <iostream>
#include <math.h>
#include <random>
#include <string>

#include <boost/math/distributions/normal.hpp>

namespace opt = boost::program_options;

void parseComLine(int argc, char **argv, std::string &infile, std::string &popmap, std::string &abcd, int &vsize, bool &two, bool &three, bool &four, int &bootstrap);
double average(std::vector<double> &v);
double stdev(std::vector<double> &v, double avg);
double calcZ(double f, double sd);

int main(int argc, char** argv) {

	std::string infile;
	std::string popmap;
	std::string abcd;
	int vsize;
	bool two = false;
	bool three = false;
	bool four = false;
	int bootstrap = 0;

	parseComLine(argc,argv,infile,popmap,abcd,vsize,two,three,four,bootstrap);

	//start random number generator
	std::default_random_engine generator;
	unsigned int seed = time(0);
	generator.seed(seed);

	fnFiles f(infile, popmap, abcd, vsize);
	f.readfiles(vsize);
	std::vector<double> returnv;
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
	double f2 = 0.0;
	double f3 = 0.0;
	double f4 = 0.0;

	if(two == true)
	{
		fs.calcF2(f,returnv);
		f2 = returnv[0];
	}
	else if (three == true)
	{
		fs.calcF3(f,returnv);
		f2 = returnv[0];
		f3 = returnv[1];
	}
	else if(four == true)
	{
		fs.calcF4(f,returnv);
		f2 = returnv[0];
		f3 = returnv[1];
		f4 = returnv[2];
	}
	else
	{
		std::cerr << "No F statistic option was selected." << std::endl;
		exit(EXIT_FAILURE);
	}
	// do bootstrapping
	std::vector<double> f2bv;
	std::vector<double> f3bv;
	std::vector<double> f4bv;

	f2bv.resize(bootstrap);
	f3bv.resize(bootstrap);
	f4bv.resize(bootstrap);

	for(int i=0; i<bootstrap; i++)
	{
		std::vector<double> bootrv;
		
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
			fsboot.calcF2(fboot,bootrv);
			f2bv[i] = bootrv[0];
		}
		else if (three == true)
		{
			fnStats fsboot(fs,bootv,f.getLength(),f.ABCDmap);
			fnFiles fboot(f,bootv,f.ABCDmap);
			fsboot.calcF3(fboot,bootrv);
			f2bv[i] = bootrv[0];
			f3bv[i] = bootrv[1];
		}
		else if(four == true)
		{
			fnStats fsboot(fs,bootv,f.getLength(),f.ABCDmap);
			fnFiles fboot(f,bootv,f.ABCDmap);
			fsboot.calcF4(fboot,bootrv);
			f2bv[i] = bootrv[0];
			f3bv[i] = bootrv[1];
			f4bv[i] = bootrv[2];
		}
		else
		{
			std::cerr << "No F statistic option was selected." << std::endl;
			exit(EXIT_FAILURE);
		}

	}
	
	double f2mean = 0.0;
	double f3mean = 0.0;
	double f4mean = 0.0;
	double f2sd = 0.0;
	double f3sd = 0.0;
	double f4sd = 0.0;

	double f2Z = 0.0;
	double f3Z = 0.0;
	double f4Z = 0.0;
	double f2Zp = 0.0;
	double f3Zp = 0.0;
	double f4Zp = 0.0;

	boost::math::normal_distribution<> zdist(0.0,1.0);
	
	if(two == true)
	{
		f2mean = average(f2bv);
		f2sd = stdev(f2bv,f2mean);
		f2Z = calcZ(f2,f2sd);
		f2Zp = 2.0*(1-boost::math::cdf(zdist,fabs(f2Z)));
		std::cout << "F2 bootstrap mean = " << f2mean << std::endl;
		std::cout << "F2 bootstrap stdev = " << f2sd << std::endl;
		std::cout << "F2 Z = " << f2Z << std::endl;
		std::cout << "F2 Z_p-val = " << f2Zp << std::endl;
	}
	else if(three == true)
	{
		f2mean = average(f2bv);
		f3mean = average(f3bv);
		f2sd = stdev(f2bv,f2mean);
		f3sd = stdev(f3bv,f3mean);
		f2Z = calcZ(f2,f2sd);
		f3Z = calcZ(f3,f3sd);
		f2Zp = 2.0*(1-boost::math::cdf(zdist,fabs(f2Z)));
		f3Zp = 2.0*(1-boost::math::cdf(zdist,fabs(f3Z)));
		std::cout << "F2 bootstrap mean = " << f2mean << std::endl;
		std::cout << "F2 bootstrap stdev = " << f2sd << std::endl;
		std::cout << "F2 Z = " << f2Z << std::endl;
		std::cout << "F2 Z_p-val = " << f2Zp << std::endl;
		std::cout << "F3 bootstrap mean = " << f3mean << std::endl;
		std::cout << "F3 bootstrap stdev = " << f3sd << std::endl;
		std::cout << "F3 Z = " << f3Z << std::endl;
		std::cout << "F3 Z_p-val = " << f3Zp << std::endl;
	}
	else if(four == true)
	{
		f2mean = average(f2bv);
		f3mean = average(f3bv);
		f4mean = average(f4bv);
		f2sd = stdev(f2bv,f2mean);
		f3sd = stdev(f3bv,f3mean);
		f4sd = stdev(f4bv,f4mean);
		f2Z = calcZ(f2,f2sd);
		f3Z = calcZ(f3,f3sd);
		f4Z = calcZ(f4,f4sd);
		f2Zp = 2.0*(1-boost::math::cdf(zdist,fabs(f2Z)));
		f3Zp = 2.0*(1-boost::math::cdf(zdist,fabs(f3Z)));
		f4Zp = 2.0*(1-boost::math::cdf(zdist,fabs(f4Z)));
		std::cout << "F2 bootstrap mean = " << f2mean << std::endl;
		std::cout << "F2 bootstrap stdev = " << f2sd << std::endl;
		std::cout << "F2 Z = " << f2Z << std::endl;
		std::cout << "F2 Z_p-val = " << f2Zp << std::endl;
		std::cout << "F3 bootstrap mean = " << f3mean << std::endl;
		std::cout << "F3 bootstrap stdev = " << f3sd << std::endl;
		std::cout << "F3 Z = " << f3Z << std::endl;
		std::cout << "F3 Z_p-val = " << f3Zp << std::endl;
		std::cout << "F4 bootstrap mean = " << f4mean << std::endl;
		std::cout << "F4 bootstrap stdev = " << f4sd << std::endl;
		std::cout << "F4 Z = " << f4Z << std::endl;
		std::cout << "F4 Z_p-val = " << f4Zp << std::endl;
	}
	else
	{
		std::cerr << "No F statistic option was selected." << std::endl;
		exit(EXIT_FAILURE);
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

double average(std::vector<double> &v)
{
	double sum = 0.0;
	for(unsigned int i=0; i<v.size(); i++)
	{
		sum+=v.at(i);
	}
	double avg = sum/(double)v.size();

	return avg;
}

double stdev(std::vector<double> &v, double avg)
{
	double var = 0.0;
	for(unsigned int i=0; i<v.size(); i++)
	{
		double dev = pow((v.at(i) - avg), 2.0);
		var+= dev;
	}
	var = var/(double)v.size();

	double sd = sqrt(var);

	return sd;
}

double calcZ(double f, double sd)
{
	double Z;
	if(sd==0)
	{
		Z=0.0;
	}
	else
	{
		Z = (0.0-f) / sd;
	}
	return Z;
}
