#ifndef POPZDSTAT_H
#define	POPZDSTAT_H

#include "Dstat.h"
#include "DstatParent.h"

#include "popZParent.h"

#include <vector>

class popZDstat: public popZParent {
public:
    popZDstat();
    void add(DstatParent *d) override;
    void calcStats() override;
    void dstatZ();
private:
    std::vector<double> D;
};

#endif	/* POPZDSTAT_H */