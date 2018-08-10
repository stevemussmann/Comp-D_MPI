/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Stats.cpp
 * Author: mussmann
 * 
 * Created on May 9, 2016, 8:03 AM
 */

#include "Stats.h"
#include <math.h>

Stats::Stats() 
{
    
}

Stats::Stats(const Stats& orig) 
{
    
}

Stats::~Stats() 
{
    
}

double Stats::chisqr(int ABBA, int BABA)
{
    int n = ABBA + BABA;
    double Xsqr = 0.0;
    if(n != 0)
    {
        double exp = (double)n/2;
    
        Xsqr+=pow((ABBA-exp), 2.0)/exp;
        Xsqr+=pow((BABA-exp), 2.0)/exp; 
    }

    return Xsqr;
}

double Stats::chisqr(double ABBA, double BABA)
{
    double n = ABBA + BABA;
    double Xsqr = 0.0;
    if(n != 0)
    {
        double exp = (double)n/2;
    
        Xsqr+=pow((ABBA-exp), 2.0)/exp;
        Xsqr+=pow((BABA-exp), 2.0)/exp; 
    }

    return Xsqr;
}

double Stats::chisqr(int L1, int L2, int L3, int L4, int R1, int R2, int R3, int R4)
{
    //uses equation for Dfoil chi-square test derived in Pease and Hahn 2015 Systematic Biology
    //it's probably more appropriate to include this in the Dfoil class, but I'm too lazy to move it.
    int L = L1+L2+L3+L4;
    int R = R1+R2+R3+R4;
    double Xsqr = 0.0;
    if(L+R != 0)
    {
        int num = pow((L-R),2);
        int denom = L+R;
        Xsqr = (double)num/(double)denom;
    }
    return Xsqr;
}

double Stats::chisqr(double L1, double L2, double L3, double L4, double R1, double R2, double R3, double R4)
{
    //uses equation for Dfoil chi-square test derived in Pease and Hahn 2015 Systematic Biology
    //it's probably more appropriate to include this in the Dfoil class, but I'm too lazy to move it.
    double L = L1+L2+L3+L4;
    double R = R1+R2+R3+R4;
    double Xsqr = 0.0;
    if(L+R != 0)
    {
        double num = pow((L-R),2);
        double denom = L+R;
        Xsqr = num/denom;
    }
    return Xsqr;
}

double Stats::average(double *arr, int length)
{
    double sum = 0.0;
    for(int i=0; i<length; i++)
    {
        sum+= arr[i];
    }
    double avg = sum/(double)length;
    
    return avg;
}

double Stats::stdev(double *arr, double avg, int length)
{
    double var = 0.0;
    for(int i=0; i<length; i++)
    {
        double dev = pow((arr[i] - avg), 2.0);
        var+= dev;
    }
    var = var/(double)length;
    
    double sd = sqrt(var);
    
    return sd;
}

double Stats::calcZ(double D, double sd)
{
    double Z;
    if(sd==0)
    {
        Z=0.0;
    }
    else
    {
        Z = (0.0 - D) / sd;
    }
    return Z;
}
