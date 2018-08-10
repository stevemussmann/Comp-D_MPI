#ifndef POPZPARTD_H
#define	POPZPARTD_H

#include "partD.h"
#include "DstatParent.h"
#include "popZParent.h"

#include <vector>

class popZpartD: public popZParent {
public:
    popZpartD();
    void add(DstatParent *d) override;
    void calcStats(std::string filename) override;
    void partdZ(std::string filename);
private:
    std::vector<double> D1;
    std::vector<double> D2;
    std::vector<double> D12;
};

#endif	/* POPZPARTD_H */
