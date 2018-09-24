# Comp-D_MPI
A program for comprehensive computation of D-statistics and population summaries using MPI. This program calculates variations of Patterson's D-statistic, including the original ABBA-BABA test as well as the partitioned-D and Dfoil tests. 

An accessory program (fn) is currently under development which will calculate the F2, F3, and F4 statistics.

## Compilation and installation

To compile, you must have a recent version of the C++ Boost Libraries and GCC 4.8.1 or newer installed on your computer. Only Boost versions 1.53 and newer have been tested. Additionally, the Boost program_options library must be present, and this is sometimes not installed automatically with the remainder of the Boost libraries. Installation of these libraries can be handled through your system's package manager to ensure all necessary components are available. For example, this is accomplished in **Ubuntu** with the command:

`sudo apt-get install libboost-dev libboost-program-options-dev libmpich2-dev autotools-dev autoconf`

Under **CentOS**, a similar command can be used:

`sudo yum install openmpi-devel boost-devel autoconf automake`

Once all dependencies are installed, download the source code, change directories into the Comp-D_MPI folder, and compile by issuing the command:

```
autoreconf
./configure
make
sudo make install
```

You can test that the program installed successfully by displaying the help menu:

`mpirun -np 1 compDmpi -h`



#### Important note for CentOS
If the configure command fails to find mpic++/mpicc, or if the  you may need to load the mpi module created by CentOS.  This is accomplished by the following command, which will need to be run in each new terminal session before executing compD:

`module load mpi`

#### Conda installation
A conda recipe for installation is under development. I am currently debugging issues encountered while compiling with MPI through conda. 

#### For Mac users
Although the conda installation does not yet work, you can still install dependencies through conda:

`conda install -c condo-forge -c anaconda -c mpi4py openmpi mpich2 autoconf automake boost boost-cpp`


## Input Files

Two files are required as input.  Examples of these files are found in the "example" folder. The first is the genotype file, which can be a phylip file of SNPs (preferred) , or a structure-formatted file.  If using a structure file, the data must be in the two-line per individual format and it should consist only of sample names with their associated genotypes.  It should not contain any of the optional data columns specifying populations, location priors, popflats, or any of the other data that can be fed to the program Structure (i.e., see the structure file output by the radseq pipeline pyRAD for an example).

The second is the taxa file, which serves as a list of the samples you wish to use for performing tests, arranged according to their position in a four- or five-taxon tree.  Sample names must match those present in the genotype file.  All samples in the taxa file must be present in the genotype file, but samples in the genotype file do not have to be included in the taxa file. See the example of the taxa file for a four-taxon test shown below:
```
48bev003 48bev004
3tsa002 3tsa010
3bcn002 3bcn003
3agr002 3agr003 3agr004 3agr001
```
In the above example, line 1 contains samples that will represent the outgroup taxon.  Line 2 represents p3 taxa,  Line 3 contains p2 taxa, and Line 4 contains p1 taxa.  One benefit Comp-D has over similar software is that it does not require you to explicitly list every combination of samples you wish to use in D-statistic calculations.  Instead, you can list each sample once on its appropriate line in the taxa file, and the software will recursively find all unique combinations of the samples.  In the example above, this results in the execution of 32 tests.

The taxa file for one of the 5-taxon tests would be formatted similarly, with lines 2 through 5 presenting p4 through p1 taxa respectively.

## Program Options for compDmpi

You can run the program to print help options with the following command:

```
mpirun -np 1 compDmpi -h
```

### Required Options
Each of these options requires input:
* **-i / --infile:** Used to specify the input file name.  See the "Acceptable File Types" heading for more information.
* **-t / --taxa:** Specify a tab delimited file to contain the taxon names upon which you will perform D tests.  See the "Taxa File Format" heading for more information.
* **-b / --bootstrap:** Specify the number of bootstrap replicates to be performed.
* **-l / --loci:** Specify the number of loci in the input file.

You must also use one of the two following boolean switches to specify the input file type:
* **-P / --phylip:** Turns on code that processes a phylip input file.
* **-s / --structure:** Turns on code that processes a structure input file.

You must also select one of the D tests to perform.  **Important: Only one of the following options can be selected each time the program is run.**  The following are the current options:
* **-d / --fourtax:** This option will perform the original four-taxon test devised by Durand et al. 2011: https://academic.oup.com/mbe/article/28/8/2239/1052492
* **-p / --partition:** This option will peform the partitioned D-statistic test developed by Eaton and Ree 2013: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3739883/
* **-f / --foil:** This option will implement the D-foil test developed by Pease and Hahn 2015: https://academic.oup.com/sysbio/article/64/4/651/1650669

For each of the D test options, the program can handle heterozygous loci in one of three ways.  The default method is to randomly select one allele whenever a heterozygous locus is encountered.  This can result in slight variations in results upon repeated runs of the program.  If you prefer a more deterministic method, one of the following two boolean switches can be used:
* **-I / --hignore:** This option will cause the program to ignore all heterozygous loci when performing D-statistic calculations.
* **-H / --hinclude:** This option will cause the program to include all heterozygote information in calculations.  This is accomplished through the SNP frequency formulas of Durand et al. 2011 and Eaton et al. 2015 (http://onlinelibrary.wiley.com/doi/10.1111/evo.12758/abstract).

### Optional commands:
* **-o / --outfile** Specify the name of the output file that contains test results (default = "outfile.txt")
* **-Z / --popstats** Specify the name of the output file that will contain population summary Z-scores (default = "popZscores.txt")

## Memory Usage:
Because this program is MPI-enabled, memory usage will increase linearly with the number of processor cores used.  Larger datasets will also use greater amounts of memory.  For example, a dataset containing 184 individuals and 179,811 SNPs was found to use ~2.9 GB per processor core. If the program crashes while calculating the first D-statistic test, it has likely run out of memory.  To verify, you can watch memory usage on your computer in real time using the `top` command to view resources being used by all processes running on your computer.

## Running the program - Example:
Imagine that you have a Phylip-formatted genotype file named "genotypes.phy" and a taxa file named "taxa.txt" which contains the taxa listed in the input file example listed above.  The file contains 14,000 SNPs, and you want to perform the four-taxon D-statistic test using the SNP frequency formulas to handle heterozygotes.  You also want to conduct 1000 bootstrap replicates for assessing significance of each test, and have 4 processor cores available to you for parallelization of the bootstrap procedure.  These options can be implemented with the following command:
```
mpirun -np 4 compDmpi -i genotypes.phy -t taxa.txt -b 1000 -l 14000 -PdH
```

## Outputs
The program writes two output files.  If you have not instructed the program to use special names for these files using the -o and -Z options, then they will be named "outfile.txt" and "popZscores.txt" respectively.

In order to ensure the output files are easily readable, you may try using the column command in Linux:
```
column -ts $'\t' outfile.txt
```

Or alternatively, open the file in Microsoft Excel or LibreOffice Calc setting the tabs as delimiters.

The first file (outfile.txt) is a tab-delimited file in which each row corresponds to a test (except for the first row, which is a header that contains column names).  The columns will vary depending upon which D-test was calculated.  For example, in the four taxon test (option -d) the first 4 columns correspond to sample names, columns 5 and 6 show the number of loci corresponding to ABBA and BABA patterns, column 7 provides the number of loci per test.  Colum 8 represents the D statistic value, 9 = standard deviation, 10 = chi-square test statistic value, 11 = chi-square p-value, 12 = Z-score test statistic value, and 13 = Z-score p-value.  Output files for Partitioned-D and D-foil tests are similar in format, but will contain several more columns corresponding to the additional site patterns and statistics that must be calculated for these tests.

As you may have already noticed, two methods are offered for assessing significance of an individual test.  Bootstrapping is necessary to calculate the standard deviation that is used to compute the Z-score (this same method is used for calculations in the pyRAD pipeline).  If you wish to avoid relying upon this method to assess statistical significance, then a chi-square test is offered as an alternative.  This option is modeled after the method Pease and Hahn 2015 use to determine significance for the D-foil test.  In practice, I find the chi-square test to more often (but not always) be a more conservative approach to assessing significance for these tests (i.e., less likely to show statistical significance).

The second file (popZscores.txt) offers a test of significance across all tests performed for a single run of compD.  This allows you to test whether a population shows significant evidence of introgression.  This computes a Z-score using the results of all tests calculated during a single run of the program, so it avoids performing the bootstrapping procedure that is necessary for the Z-score calculations you find in the outfile.txt file.  **Important: In order for this option to produce a meaningful result, you must input samples for the taxa of interest that represent members of the same populations.**

## List of bug fixes, additions, and other changes:
2018-09-24:
* Added command line option to compD which allows the user to specify the number of extra columns in a Structure-formatted file. This will allow Structure input to be more flexible.
* The program fn is still in development.

2018-08-10:
* Switched to compilation with autotools. Conda installation recipe is under development

2018-07-27:
* Began adding accessory program to calculate F2, F3, and F4 statistics.  Program (fn) compiles but is not yet functional.  User guide will be forthcoming.

2018-06-28:
* Major bug fix - p-values associated with Z-scores are now correct.  If you were relying upon the Z-scores for determining statistical significance, then you should re-run any analyses done prior to this date.  Under certain conditions the Z-score was being converted from a double to integer prior to calculation of the p-value, leading to incorrect p-values.  The p-values associated with the Chi-squared test were unaffected by this bug.

## Legacy Code for compDmpi
I originally wrote the program with the intention of using the whole sequence of a RAD locus for consideration when comparing alleles.  However, I found this often produced messy results as well as increasing computation time and memory usage.  Therefore, I consider the remaining options to be "deprecated" and generally advise against their use.  However, this functionality is still included in the program if you want to experiment with these options.  

To implement these methods, do not specify a file format option (-P or -s).  This will make the program you are inputting sequences in a pyRAD .alleles file format. I also experimented with a few options that relate to whole loci, and these options are listed below.  If you want to ignore certain sites within sequences (such as gaps or uncalled bases) you can use the following options:
* **-n / --nremove** Boolean switch.  When used, it will force the program to ignore uncalled bases (N) within alleles when performing string matching operations.
* **-g / --gap** Boolean switch.  When used, it will force the program to ignore any gaps (-) within alleles when performing string matching operations.
