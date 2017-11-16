#ifndef POPZPARENT_H
#define	POPZPARENT_H

#include "DstatParent.h"

#include <vector>

class popZParent: public Stats {
public:
    popZParent();
    //popZParent(const popZParent& orig);
    //virtual ~popZParent();
    virtual void add(DstatParent *d) = 0;
    virtual void calcStats() = 0;
protected:
    double* toArr(std::vector<double> &vec, unsigned int length);
};

#endif	/* POPZPARENT_H */