/* 
 * File:   locusfile.h
 * Author: Steve
 *
 * Created on July 7, 2015, 9:03 PM
 */

#ifndef LOCUSFILE_H
#define	LOCUSFILE_H

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

class locusfile{
public:
    locusfile(int vectorsize);
    void readInput(std::string infile, int locnumber);
    void readInput(std::string infile, int locnumber, int offset);
    void readInput(std::string infile, int locnumber, bool phylip);
    void removeN(int locnumber, bool Nremoveflag, bool gapignoreflag);
    std::vector<int> findInformative(locusfile &current, int locnumber, std::unordered_map <std::string,int> indlist, int ntaxa, bool hetIgnore);
    std::string iupac(std::string ambig);
    void AddData(int i, std::string name, std::string sequence, int allele);
    void AddData(int i, std::string name, std::string allele1, std::string allele2, std::vector<std::unordered_map<std::string,int> > freq);
    int GetSize(int i);
    int GetSeqSize(int i);
    void calcFreq(int locnumber, std::vector<std::vector<std::string>> &taxa);
    std::string getMajorAllele(int locus, int taxon, int my_rank);
    std::string GetSeq(int i, int j);
    std::string GetName(int i, int j);
    virtual ~locusfile();
private:
    int vectorsize;
    std::vector<std::vector<std::string> > species;
    std::vector<std::vector<std::string> > alleles;
    std::vector<std::vector<std::unordered_map<std::string,int> > > freqs;
    std::vector<int> findLocation(std::string sequence, char findN, char findDash, bool Nremoveflag, bool gapignoreflag);
};

#endif	/* LOCUSFILE_H */

