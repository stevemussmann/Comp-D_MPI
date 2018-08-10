#ifndef POPZDFOIL_H
#define	POPZDFOIL_H

#include "Dfoil.h"
#include "DstatParent.h"

#include "popZParent.h"

#include <vector>

class popZDfoil: public popZParent {
public:
    popZDfoil();
    void add(DstatParent *d) override;
    void calcStats(std::string filename) override;
    void dfoilZ(std::string filename);
private:
    std::vector<double> DFO;
    std::vector<double> DIL;
    std::vector<double> DFI;
    std::vector<double> DOL;
};

#endif	/* POPZDFOIL_H */
