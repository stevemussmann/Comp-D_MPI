#include "popZParent.h"
#include "DstatParent.h"

popZParent::popZParent(){

}

double* popZParent::toArr(std::vector<double> &vec, unsigned int length){
    double *arr = new double[length];
    for(unsigned int i=0; i<length; i++){
	arr[i] = vec.at(i);
    }
    
    return arr;
}

