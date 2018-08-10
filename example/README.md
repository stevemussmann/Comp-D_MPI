# Example Files

## Genotype Files
* example.phy = Phylip-formatted file containing 11899 SNPs from 112 ingroup samples representing 6 genetic populations, and 5 outgroup samples
* example.str = Structure-formatted file containing the same data as example.phy

## Taxa File
* taxafile.txt = example taxa file for the four-taxon D-statistic test

## Example files for usage in fn
* popmap.txt = tab-delimited population map that stores the population information for the example.phy and example.str files.
* abcd_map.txt = tab-delimited file for calculating the F4, F3, and F2 statistics in fn
* abc_map.txt = tab-delimited file for calculating the F3 and F2 statistics in fn
* ab_map.txt = tab-delimited file for calculating the F2 statistic in fn

## Example commands
To run the four-taxon test on the contents of taxafile.txt using the Phylip-formatted file:

`mpirun -np 4 compDmpi -i example.phy -t taxafile.txt -l 11899 -b 100 -o test.out -Z test.popZ.out -Pd`

To run the same command on taxafile.txt using the Structure-formatted file:

`mpirun -np 4 compDmpi -i example.str -t taxafile.txt -l 11899 -b 100 -o test.out -Z test.popZ.out -sd`
