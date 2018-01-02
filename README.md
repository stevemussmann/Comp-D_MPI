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

List of current required options:
* **-i / --infile:** Used to specify the input file name.  See the "Acceptable File Types" heading for more information.
* **-t / --taxa:** Specify a tab delimited file to contain the taxon names upon which you will perform D tests.  See the "Taxa File Format" heading for more information.
* **-b / --bootstrap:** Specify the number of bootstrap replicates to be performed.
* **-l / --loci:** Specify the number of loci in the input file.  

You must also select one of the D tests to perform.  Important: Only one of the following options can be selected each time the program is run.  The following are the current options:
* **-d / --fourtax:** This option will perform the original four-taxon test devised by Durand et al. 2011: https://academic.oup.com/mbe/article/28/8/2239/1052492
* **-p / --partition:** This option will peform the partitioned D-statistic test developed by Eaton and Ree 2013: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3739883/

