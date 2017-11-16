/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Stats.h
 * Author: mussmann
 *
 * Created on May 9, 2016, 8:03 AM
 */

#ifndef STATS_H
#define STATS_H

class Stats {
    public:
	Stats();
	Stats(const Stats& orig);
	virtual ~Stats();
    protected:
        double chisqr(int ABBA, int BABA);
	double chisqr(double ABBA, double BABA);
        double chisqr(int L1, int L2, int L3, int L4, int R1, int R2, int R3, int R4);
	double chisqr(double L1, double L2, double L3, double L4, double R1, double R2, double R3, double R4);
	double calcZ(double D, double sd);
	double average(double *arr, int length);
	double stdev(double *arr, double avg, int length);
};

#endif /* STATS_H */