/* 
 * File:   locus.cpp
 * Author: Steve
 * 
 * Created on July 15, 2015, 7:38 AM
 */

#include "locus.h"
#include <string>
#include <vector>

locus::locus()
{
    species;
    sequences;
}

void locus::AddData(std::string name, std::string sequence, int allele)
{
    if(allele == 0)
    {
        species.push_back(name);
    }
    //put the sequence onto the vector
    sequences.push_back(sequence);
}

void locus::AddData(std::string name, std::string allele1, std::string allele2)
{
    species.push_back(name);
    sequences.push_back(allele1);
    sequences.push_back(allele2);
}

int locus::GetSize()
{
    int length = species.size();
    
    return length;
}

int locus::GetSeqSize()
{
    int length = sequences.size();
    
    return length;
}

std::string locus::GetName(int i)
{
    std::string name = species.at(i);
    
    return name;
}

std::string locus::GetSeq(int i)
{
    std::string seq = sequences.at(i);
    
    return seq;
}

void locus::SetSeq(std::string seq, int i)
{
    sequences.at(i) = seq;
}

locus::locus(const locus& orig) {
}

locus::~locus() {
}

