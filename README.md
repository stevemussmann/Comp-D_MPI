# Comp-D_MPI
A program for comprehensive computation of D-statistics and population summaries using MPI

This program calculates variations of Patterson's D-statistic, including the original ABBA-BABA test as well as the partitioned-D and Dfoil tests.  A full user guide for the program is forthcoming.

## Compilation

To compile, you must have a recent version of the C++ Boost Libraries and GCC 4.8.1 or newer installed on your computer. Only Boost versions 1.54 and newer have been tested. Download the source code, change directories into the Comp-D_MPI folder, and compile by issuing the command:

`make`

If you wish to automatically install the program to your path, then issue the command:

`sudo make install`

This automatically installs the program into the /usr/local/bin directory.

## Program Options

You can run the program to print help options with the following command:

```
compDmpi -h
```

###Required Options
List of current required options:
* **-i / --infile:** Used to specify the input file name.  See the "Acceptable File Types" heading for more information.
* **-t / --taxa:** Specify a tab delimited file to contain the taxon names upon which you will perform D tests.  See the "Taxa File Format" heading for more information.
* **-b / --bootstrap:** Specify the number of bootstrap replicates to be performed.
* **-l / --loci:** Specify the number of loci in the input file.

If you are inputting a structure or phylip file, you must use one of the two following boolean switches:
* **-P / --phylip:** Turns on code that processes a phylip input file.  Given that this file format is more standardized and less flexible than the structure format, I consider this to be the "preferred" file format choice for input.
* **-s / --structure:** Turns on code that processes a structure input file.  This is the two-line per individual format.  The data in the file should consist only of genotypes.  It should not contain any extra columns specifying populations, popflag, or any of the other extra  data that can be fed to the program Structure (i.e., see the structure file output by the radseq pipeline pyRAD for an example).

You must also select one of the D tests to perform.  Important: Only one of the following options can be selected each time the program is run.  The following are the current options:
* **-d / --fourtax:** This option will perform the original four-taxon test devised by Durand et al. 2011: https://academic.oup.com/mbe/article/28/8/2239/1052492
* **-p / --partition:** This option will peform the partitioned D-statistic test developed by Eaton and Ree 2013: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3739883/
* **-f / --foil:** This option will implement the D-foil test developed by Pease and Hahn 2015: https://academic.oup.com/sysbio/article/64/4/651/1650669

For each of these methods, the program can handle heterozygous loci in one of three ways.  The default method is to randomly select one allele whenever a heterozygous locus is encountered.  This can result in slight variations in results upon repeated runs of the program.  If you prefer a more deterministic method, one of the following two boolean switches can be used:
* **-I / --hignore:** This option will cause the program to ignore all heterozygous loci when performing D-statistic calculations.
* **-H / --hinclude:** This option will cause the program to include all heterozygote information in calculations.  

Optional commands:
* **-o / --outfile** Specify the name of the output file that contains test results (default = "outfile.txt")



## Legacy Code
I originally wrote the program with the intention of using the whole sequence of a RAD locus for consideration when comparing alleles.  However, I found this often produced messy results as well as increasing computation time and memory usage.  Therefore, I consider the remaining options to be "deprecated" and generally advise against their use.  However, this functionality is still included in the program if you want to experiment with these options.  

To implement these methods, do not specify a file format option (-P or -s).  This will make the program you are inputting sequences in a pyRAD .alleles file format. I also experimented with a few options that relate to whole loci, and these options are listed below.  If you want to ignore certain sites within sequences (such as gaps or uncalled bases) you can use the following options:
* **-n / --nremove** Boolean switch.  When used, it will force the program to ignore uncalled bases (N) within alleles when performing string matching operations.
* **-g / --gap** Boolean switch.  When used, it will force the program to ignore any gaps (-) within alleles when performing string matching operations.
