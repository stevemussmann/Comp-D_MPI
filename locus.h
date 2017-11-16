/* 
 * File:   locus.h
 * Author: Steve
 *
 * Created on July 15, 2015, 7:38 AM
 */

#ifndef LOCUS_H
#define	LOCUS_H

#include <string>
#include <vector>

class locus {
public:
    locus();
    void AddData(std::string name, std::string sequence, int allele);
    void AddData(std::string name, std::string allele1, std::string allele2);
    int GetSize();
    int GetSeqSize();
    std::string GetName(int i);
    std::string GetSeq(int i);
    void SetSeq(std::string seq, int i);
    locus(const locus& orig);
    virtual ~locus();
protected:
    std::vector<std::string> species;
    std::vector<std::string> sequences;
};

#endif	/* LOCUS_H */

