# Comp-D_MPI
A program for comprehensive computation of D-statistics and population summaries using MPI

This program calculates variations of Patterson's D-statistic, including the original ABBA-BABA test as well as the partitioned-D and Dfoil tests.  A full user guide for the program is forthcoming.

## Compilation

To compile, you must have a recent version of the C++ Boost Libraries and GCC 4.8.1 or newer installed on your computer. Only Boost versions 1.54 and newer have been tested. Download the source code, change directories into the Comp-D_MPI folder, and compile by issuing the command:

`make`

If you wish to automatically install the program to your path, then issue the command:

`sudo make install`

This automatically installs the program into the /usr/local/bin directory.

## Input Files

Two files are required as input.  The first is the genotype file, which can be a phylip file of SNPs (preferred) , or a structure-formatted file.  If using a structure file, the data must be in the two-line per individual format and it should consist only of sample names with their associated genotypes.  It should not contain any of the optional data columns specifying populations, location priors, popflats, or any of the other data that can be fed to the program Structure (i.e., see the structure file output by the radseq pipeline pyRAD for an example).

The second is the taxa file, which serves as a list of the samples you wish to use for performing tests, arranged according to their position in a four- or five-taxon tree.  Sample names must match those present in the genotype file.  All samples in the taxa file must be present in the genotype file, but samples in the genotype file do not have to be included in the taxa file. See the example of the taxa file for a four-taxon test shown below:
```
48bev003 48bev004
3tsa002 3tsa010
3bcn002 3bcn003
3agr002 3agr003 3agr004 3agr001
```
In the above example, line 1 contains samples that will represent the outgroup taxon.  Line 2 represents p3 taxa,  Line 3 contains p2 taxa, and Line 4 contains p1 taxa.  One benefit Comp-D has over similar software is that it does not require you to explicitly list every combination of samples you wish to use in D-statistic calculations.  Instead, you can list each sample once on its appropriate line in the taxa file, and the software will recursively find all unique combinations of the samples.  In the example above, this results in the execution of 32 tests.

The taxa file for one of the 5-taxon tests would be formatted similarly, with lines 2 through 5 presenting p4 through p1 taxa respectively.

## Program Options

You can run the program to print help options with the following command:

```
compDmpi -h
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

## Running the program - Example:

Imagine that you have a Phylip-formatted genotype file named "genotypes.phy" and a taxa file named "taxa.txt" which contains the taxa listed in the input file example listed above.  The file contains 14,000 SNPs, and you want to perform the four-taxon D-statistic test using the SNP frequency formulas to handle heterozygotes.  You also want to conduct 1000 bootstrap replicates for assessing significance of each test, and have 4 processor cores available to you for parallelization of the bootstrap procedure.  These options can be implemented with the following command:
```
mpirun -np 4 compDmpi -i genotypes.phy -t taxa.txt -b 1000 -l 14000 -PdH
```

## Legacy Code
I originally wrote the program with the intention of using the whole sequence of a RAD locus for consideration when comparing alleles.  However, I found this often produced messy results as well as increasing computation time and memory usage.  Therefore, I consider the remaining options to be "deprecated" and generally advise against their use.  However, this functionality is still included in the program if you want to experiment with these options.  

To implement these methods, do not specify a file format option (-P or -s).  This will make the program you are inputting sequences in a pyRAD .alleles file format. I also experimented with a few options that relate to whole loci, and these options are listed below.  If you want to ignore certain sites within sequences (such as gaps or uncalled bases) you can use the following options:
* **-n / --nremove** Boolean switch.  When used, it will force the program to ignore uncalled bases (N) within alleles when performing string matching operations.
* **-g / --gap** Boolean switch.  When used, it will force the program to ignore any gaps (-) within alleles when performing string matching operations.
