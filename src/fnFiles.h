#ifndef FNFILES_H
#define FNFILES_H

#include <string>
#include <unordered_map>
#include <vector>

class fnFiles {
	public:
		fnFiles(std::string i, std::string p, std::string a, int vectorsize);
		fnFiles(fnFiles f, std::vector<int> &v, std::unordered_map<std::string,std::string> m);
		void readfiles(int vectorsize, bool p, bool s, std::string missing, int offset);
		unsigned int getLength();
		std::unordered_map <std::string,int> getLocus(std::string s, int i);
		std::unordered_map <std::string,std::string> ABCDmap;
		void checkF4();
		void checkF3();
		void checkF2();
	private:
		std::string infile;
		std::string popfile;
		std::string outgroup;
		std::string ABCDfile;
		void readPhylip();
		void readStructure(int offset, int l, std::string m);
		void readPopfile();
		void readABCDfile();
		void blacklist();
		std::unordered_map<std::string,std::vector<std::unordered_map<std::string,int > > > data;
		std::unordered_map <std::string,std::string> popmap;
		std::string iupac(std::string ambig);
		void checkPops(std::vector<std::string> &v);
};

#endif
