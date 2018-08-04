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
		std::unordered_map <std::string,int> getDLocus(int i);
		std::unordered_map <std::string,int> getALocus(int i);
		std::unordered_map <std::string,int> getBLocus(int i);
		std::unordered_map <std::string,int> getCLocus(int i);
		std::unordered_map <std::string,int> getLocus(std::string s, int i);
		std::unordered_map <std::string,std::string> ABCDmap;
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
		std::string iupac(std::string ambig);
};

#endif
