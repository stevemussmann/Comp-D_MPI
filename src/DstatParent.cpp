/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   DstatParent.cpp
 * Author: Steve
 * 
 * Created on May 9, 2016, 8:03 AM
 */

#include "DstatParent.h"
#include "Stats.h"
#include "locusfile.h"

#include <math.h>
#include <fstream>
#include <iostream>

DstatParent::DstatParent() 
{
    
}

DstatParent::DstatParent(const DstatParent& orig) 
{
    
}

DstatParent::~DstatParent() 
{
    
}

void DstatParent::writeout(std::string *array, std::string output, int i, bool hetIgnore, bool hetInclude )
{
    std::ofstream outfile;
    if(i == 0)
    {
        outfile.open(output, std::ios::out );
    }
    else
    {
        outfile.open(output, std::ios::out | std::ios::app );
    }
    
    if(outfile.is_open())
    {
	this->write(array, outfile, i, hetIgnore, hetInclude);
    }
    else
    {
        std::cerr << "Unable to open file for output" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    outfile.close();
}
