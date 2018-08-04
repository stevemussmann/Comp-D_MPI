#ifndef FNFILES_H
#define FNFILES_H

#include <string>
#include <unordered_map>
#include <vector>

class fnFiles {
	public:
		//fnFiles(std::string i, std::string p, std::string o, std::string a, int vectorsize);
		fnFiles(std::string i, std::string p, std::string a, int vectorsize);
		void readfiles();
		unsigned int getLength();
		std::unordered_map <std::string,int> getOutgroupLocus(int i);
		std::unordered_map <std::string,int> getALocus(int i);
		std::unordered_map <std::string,int> getBLocus(int i);
		std::unordered_map <std::string,int> getCLocus(int i);
	private:
		std::string infile;
		std::string popfile;
		std::string outgroup;
		std::string ABCDfile;
		void readPhylip();
		void readPopfile();
		void readABCDfile();
		void blacklist();
		std::unordered_map<std::string,std::vector<std::unordered_map<std::string,int > > > data;
		std::unordered_map <std::string,std::string> popmap;
		std::unordered_map <std::string,std::string> ABCDmap;
		std::vector<std::unordered_map <std::string,int> > A;
		std::vector<std::unordered_map <std::string,int> > B;
		std::vector<std::unordered_map <std::string,int> > C;
		std::vector<std::unordered_map <std::string,int> > D;
		std::string iupac(std::string ambig);
};

#endif
