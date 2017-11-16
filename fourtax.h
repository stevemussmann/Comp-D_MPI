/* 
 * File:   fourtax.h
 * Author: Steve
 *
 * Created on July 7, 2015, 9:03 PM
 */

#ifndef FOURTAX_H
#define	FOURTAX_H

#include "locusfile.h"
#include <string>
#include <vector>
#include <unordered_map>

class fourtax {
public:
    fourtax(int locusnum, int taxnum);
    fourtax(const fourtax& orig);
    void calculatePattern(int locus, int ntaxa, locusfile &file, int my_rank);
    std::string getPattern(int i);
    double getFreq(int locus, int taxon);
    void populateDtest(std::vector<int> &keep, locusfile &file, std::unordered_map <std::string,int> &indlist, std::default_random_engine &generator, int my_rank, int ntaxa);
    void populateDtest(std::vector<int> &keep, locusfile &file, std::unordered_map<std::string,int> &indlist, int my_rank, int ntaxa);
    virtual ~fourtax();
private:
    std::vector<int> locusnum;
    std::vector<int> taxnum;
    std::vector<int> number;
    std::vector<std::vector<std::string> > major;
    std::vector<std::vector<std::string> > names;
    std::vector<std::vector<std::string> > seq;
    std::vector<std::vector<double> > freq;
    std::vector<std::string> pattern;
    
};

#endif	/* FOURTAX_H */

