/* 
 * File:   dstat_main.cpp
 * Author: Steve
 *
 * Created on July 7, 2015, 8:00 PM
 */

#include <algorithm> //std::unique, std::distance
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string> //std::string
#include <unordered_map> //std::unordered_map
#include <vector> //std::vector

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/program_options.hpp>

#include "mpi.h"

#include "locusfile.h"
#include "fourtax.h"
#include "Dstat.h"
#include "partD.h"
#include "Dfoil.h"
#include "DstatParent.h"

#include "popZParent.h"
#include "popZDstat.h"
#include "popZDfoil.h"
#include "popZpartD.h"

using namespace std;

namespace opt = boost::program_options;

void readTaxa(string infile, std::vector<vector<string>> &taxa);
void combinations(vector<vector<string> > &array, unsigned int i, 
        vector<string> accum, vector<vector<string> > &comb);
void parseComLine(string &infile, string &taxafile, int &locnumber, 
        int &bootstrap, string &output, int argc, char **argv, 
        bool &fourtaxflag, bool &partDflag, bool &Dfoilflag, 
        bool &Nremoveflag, bool &gapignoreflag, bool &strflag, bool &phylip, 
	bool &hetIgnore, bool &hetInclude, string &popoutput, int &my_rank);


int main(int argc, char** argv) {
    
    int my_rank; //initialize processor rank for MPI
    int p; //initialize number of processors
    
    MPI_Init(&argc, &argv); //start MPI
    
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); //assign values to my_rank
    
    MPI_Comm_size(MPI_COMM_WORLD, &p); //assign value to 
    
    double start = 0.0;
    if(my_rank == 0)
    {
        start = MPI_Wtime();
    }
    
    string infile;
    string taxafile;
    //int test;
    int locnumber;
    int bootstrap;
    string output;
    string popoutput;
    bool fourtaxflag = false; //turns the four taxon test on/off
    bool partDflag = false; //turns the partitioned-D test on/off
    bool Dfoilflag = false; //turns the D-FOIL test on/off
    bool Nremoveflag = false; //turns the N removal function on/off
    bool gapignoreflag = false; //turns the gap removal function on/off
    bool hetIgnore = false; //turns off any use of heterozygotes instead of randomly picking an allele
    bool hetInclude = false; //uses all heterozygote information for calculations
    bool strflag = false; //turn on to input a structure file
    bool phylip = false; //turn on to input a phylip file
    
    //parse the command line
    parseComLine(infile, taxafile, locnumber, bootstrap, output, argc, argv, 
		 fourtaxflag, partDflag, Dfoilflag, Nremoveflag, gapignoreflag, 
		 strflag, phylip, hetIgnore, hetInclude, popoutput, my_rank);
    
    //calculate number of bootstraps per processor
    int mpiboot = bootstrap/p;
    int remainder = bootstrap%p;
    if(remainder>0)
    {
        mpiboot+=1;
    }
    bootstrap = mpiboot*p;
    
    //start random number generator
    default_random_engine generator;
    unsigned int seed = (time(0) + my_rank);
    generator.seed(seed);
    
    vector<vector<string>> taxa;
    
    clock_t readtax;
    readtax = clock();
    //read input taxa list file
    readTaxa(taxafile, taxa);
    double readtaxrun = (clock() - readtax) / (double)CLOCKS_PER_SEC;
    cout << "Time to read taxa file on proc " << my_rank << " = " << readtaxrun << " seconds." << endl;
    
    MPI_Barrier(MPI_COMM_WORLD); //make sure all have read taxa file before proceeding
    
    locusfile newfile(locnumber); //allocate memory to hold contents of input file
    
    clock_t readalleles;
    readalleles = clock();
    
    if(strflag == true)
    {
        newfile.readInput(infile, locnumber, 0);
        
    }
    else if( phylip == true)
    {
        newfile.readInput(infile, locnumber, phylip);
    }
    else
    {
        newfile.readInput(infile, locnumber); //read the input file
    }
    
    
    double readallelesrun = (clock() - readalleles) / (double)CLOCKS_PER_SEC;
    cout << "Time to read alleles file on proc " << my_rank << " = " << readallelesrun << " seconds." << endl;
    
    MPI_Barrier(MPI_COMM_WORLD); //make sure all have read alleles file before proceeding
    
    //do error correction (i.e., remove Ns from input file)
    if(Nremoveflag==true || gapignoreflag==true)
    {
        newfile.removeN(locnumber, Nremoveflag, gapignoreflag);
    }
    
    //get all unique combinations
    vector<string> accum;
    vector<vector<string> > comb;
    unsigned int i=0;
    combinations(taxa, i, accum, comb);

    popZParent *pop;
    
    if(fourtaxflag==true)
    {
	popZDstat *popZD = new popZDstat();
	pop = popZD;
    }
    else if(partDflag==true)
    {
	popZpartD *popZD = new popZpartD();
	pop = popZD;
    }
    else if(Dfoilflag==true)
    {
	popZDfoil *popZD = new popZDfoil();
	pop = popZD;
    }
    else
    {
	cout << "Invalid test type." << endl;
	exit(EXIT_FAILURE);
    }
    
    //get lists of p3 taxa and outgroups
    newfile.calcFreq(locnumber, taxa);
    
    
    //score alleles for loci
    
    
    
    //cout << "test type declared." << endl << endl;
	for(vector<vector<string>>::size_type i=0; i<comb.size(); i++)
	{
		DstatParent *test;
		int ntaxa;
		if(fourtaxflag==true)
		{
			Dstat *Dstattest = new Dstat();
			test = Dstattest;
			ntaxa = 4;
		}
		else if(partDflag==true)
		{
			partD *partDtest = new partD();
			test = partDtest;
			ntaxa = 5;
		}
		else if(Dfoilflag==true)
		{
			Dfoil *Dfoiltest = new Dfoil();
			test = Dfoiltest;
			ntaxa = 5;
		}
		else
		{
			cout << "Invalid test type." << endl;
			exit(EXIT_FAILURE);
		}
		
		unordered_map <string,int> indlist; //hash of individuals
		string *indarray; //list of individuals for writing output
		indarray = new string[comb[i].size()];
		
		for(vector<string>::size_type j=0; j<comb[i].size(); j++ ) //construct hash of individuals
		{
			indlist[comb[i][j]] = j;
			indarray[j] = comb[i][j];
		}
		
		locusfile current(locnumber); //declare new locusfile to hold data
		vector<int> keep = newfile.findInformative(current, locnumber, indlist, ntaxa, hetIgnore);
		
		//Dstat Dstattest;
		if(my_rank == 0)
		{
			fourtax dtest(keep.size(), ntaxa);
			if(hetIgnore==true || (hetIgnore==false && hetInclude==false))
			{
				dtest.populateDtest(keep, current, indlist, generator, my_rank, ntaxa);
				test->calcDs(dtest, keep.size(), ntaxa, current, my_rank);
				test->calcStats(keep.size());
			}
			else if(hetInclude==true)
			{
				dtest.populateDtest(keep, current, indlist, my_rank, ntaxa);
				test->calcPolyDs(dtest, keep.size());
				test->polyCalcStats(keep.size());
			}
			else
			{
				cerr << "This code should not be reachable." << endl;
				exit(EXIT_FAILURE);
			}
		}
		//calculate proportion of discordant loci
		MPI_Barrier(MPI_COMM_WORLD);
		
		test->bootstrap(mpiboot, bootstrap, indlist, generator, ntaxa, current, keep, my_rank, i, comb.size(), indarray, output, hetIgnore, hetInclude);
		
		if(my_rank == 0){
			//put D stats into arrays
			pop->add(test);
		}
		
		delete test;
	}

    MPI_Barrier(MPI_COMM_WORLD);
    
    if(my_rank == 0)
    {
        double stop = MPI_Wtime();
        double runtime = stop-start;
        cout << "Time to completion was " << runtime << " seconds." << endl;
    }
    
    if(my_rank == 0)
    {
	pop->calcStats(popoutput);
	std::cout << "Population summary Z-scores were written to " << popoutput << std::endl;
    }
    
    delete pop;
    
    MPI_Finalize();
    
    return 0;
}

void readTaxa(string infile, vector<vector<string>> &taxa)
{    
    ifstream taxfile(infile);
    if( taxfile.is_open())
    {
        string line;
        while(std::getline(taxfile, line))
        {
            if(line.empty())
            {
                cout << "Empty line in taxa file. Check input file if program fails" << endl;
            }
            else
            {
                stringstream ss(line);
                string buf;
                vector<string> vec;
                while(ss >> buf)
                {
                    vec.push_back(buf);
                }   
            taxa.push_back(vec);
            }
        }
    }
    else
    {
        cerr << "Unable to open " << infile << endl;
        exit(EXIT_FAILURE);
    }
}

void combinations(vector<vector<string> > &array, unsigned int i, vector<string> accum, vector<vector<string> > &comb)
{
    if (i == array.size())
    {
        comb.push_back(accum);
    }
    else
    {
        vector<string> row = array[i];
        for(unsigned int j = 0; j < row.size(); j++)
        {
            vector<string> tmp(accum);
            tmp.push_back(row[j]);
            combinations(array,i+1,tmp,comb);
        }
    }
}

void parseComLine(string &infile, string &taxafile, int &locnumber, int &bootstrap, 
        string &output, int argc, char **argv, bool &fourtaxflag, bool &partDflag, 
        bool &Dfoilflag, bool &Nremoveflag, bool &gapignoreflag, bool &strflag, 
        bool &phylip, bool &hetIgnore, bool &hetInclude, string &popoutput, int &my_rank)
{
    opt::options_description desc("--- Option Descriptions ---");
    desc.add_options()
            ("help,h", "Prints this help message.")
            ("infile,i", opt::value<string>(&infile)->required(), "Specifies the input alleles file name.")
            ("taxa,t", opt::value<string>(&taxafile)->required(), "Specifies the input list of taxa.")
            ("bootstrap,b", opt::value<int>(&bootstrap)->required(), "Specifies the number of bootstrap replicates to be performed.")
            ("loci,l", opt::value<int>(&locnumber)->required(), "Specifies the number of loci in the input file.")
            ("fourtax,d", opt::bool_switch(&fourtaxflag), "Turns on the 4-taxon Test.")
            ("partition,p", opt::bool_switch(&partDflag), "Turns on the Partitioned-D Test.")
            ("foil,f", opt::bool_switch(&Dfoilflag), "Turns on the Dfoil Test.")
            ("gap,g", opt::bool_switch(&gapignoreflag), "Turns on the function to ignore gaps in sequences.")
            ("nremove,n", opt::bool_switch(&Nremoveflag), "Turns on the function to remove Ns from sequences.")
	    ("hignore,I", opt::bool_switch(&hetIgnore), "Turns on function to ignore any heterozygous loci.")
	    ("hinclude,H", opt::bool_switch(&hetInclude), "Turns on function to include all heterozygote information in calculations.")
            ("structure,s", opt::bool_switch(&strflag), "Use this option to input a structure file.")
            ("phylip,P", opt::bool_switch(&phylip), "Use this option to input a phylip file.")
            ("outfile,o", opt::value<string>(&output)->default_value("outfile.txt"), "Specifies the name of the output file.")
            ("popstats,Z", opt::value<string>(&popoutput)->default_value("popZscores.txt"), "Specifies the name of the output file containing population summary Z scores.")
    ;
    
    opt::variables_map vm;
    try
    {
        opt::store(opt::parse_command_line(argc, argv, desc), vm);
    
        if(vm.count("help"))
        {
            if(my_rank == 0)
            {
                cout << "This program was created for calculating D-statistics given " << endl;
                cout << "a file of alleles produced by the ddRAD pipeline pyRAD." << endl;
                cout << desc << endl;
            }
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        
        opt::notify(vm); // throws an error if there are any problems
    }
    catch(opt::required_option& e) //catch errors resulting from required options
    {
        if(my_rank == 0)
        {
            cerr <<endl << "ERROR: " << e.what() << endl << endl;
            cout << desc << endl;
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    catch(opt::error& e) // catch other command line errors
    {
        if(my_rank == 0)
        {
            cerr << endl << "ERROR: " << e.what() << endl << endl;
            cout << desc << endl;
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    if((fourtaxflag==true && partDflag==true) || (fourtaxflag==true && Dfoilflag==true))
    {
        if(my_rank == 0)
        {
            cerr << "Cannot run four taxon and five taxon tests simultaneously." << endl;
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    if(fourtaxflag==false && partDflag==false && Dfoilflag==false)
    {
        if(my_rank == 0)
        {
            cerr << "Must select at least one test to run." << endl;
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    if(hetIgnore==true && hetInclude==true)
    {
	if(my_rank == 0)
	{
	    cerr << "Cannot use options to ingore and include heterozygotes simultaneously." << endl;
	}
	MPI_Finalize();
	exit(EXIT_FAILURE);
    }
    
}
