# Comp-D_MPI
A program for comprehensive computation of D-statistics and population summaries using MPI

This program calculates variations of Patterson's D-statistic, including the original ABBA-BABA test as well as the partitioned-D and Dfoil tests.  A full user guide for the program is forthcoming.

## Compilation

To compile, you must have a recent version of the C++ Boost Libraries and GCC 4.8.1 or newer installed on your computer. Only Boost versions 1.54 and newer have been tested. Download the source code, change directories into the Comp-D_MPI folder, and compile by issuing the command:

`make`

If you wish to automatically install the program to your path, then issue the command:

`sudo make install`

This automatically installs the program into the /usr/local/bin directory.
