all : compDmpi

Stats.o : Stats.h Stats.cpp
	g++ -O3 -std=c++11 -Wall -c Stats.cpp
	
DstatParent.o : DstatParent.h DstatParent.cpp Stats.h
	g++ -O3 -std=c++11 -Wall -c DstatParent.cpp

locusfile.o : locusfile.h locusfile.cpp
	mpic++ -O3 -std=c++11 -Wall -c locusfile.cpp

fourtax.o : fourtax.h fourtax.cpp
	g++ -O3 -std=c++11 -Wall -c fourtax.cpp

Dstat.o : Dstat.h Dstat.cpp fourtax.h DstatParent.h Stats.h
	mpic++ -O3 -std=c++11 -Wall -c Dstat.cpp

partD.o : partD.h partD.cpp fourtax.h DstatParent.h Stats.h
	mpic++ -O3 -std=c++11 -Wall -c partD.cpp

Dfoil.o : Dfoil.h Dfoil.cpp fourtax.h DstatParent.h Stats.h
	mpic++ -O3 -std=c++11 -Wall -c Dfoil.cpp

popZParent.o : popZParent.h popZParent.cpp Stats.h Dstat.h partD.h Dfoil.h
	g++ -O3 -std=c++11 -Wall -c popZParent.cpp

popZDstat.o : popZDstat.h popZDstat.cpp Stats.h Dstat.h DstatParent.h
	g++ -O3 -std=c++11 -Wall -c popZDstat.cpp

popZpartD.o : popZpartD.h popZpartD.cpp Stats.h Dstat.h DstatParent.h
	g++ -O3 -std=c++11 -Wall -c popZpartD.cpp
	
popZDfoil.o : popZDfoil.h popZDfoil.cpp Stats.h Dstat.h DstatParent.h
	g++ -O3 -std=c++11 -Wall -c popZDfoil.cpp
	
dstat_main.o : dstat_main.cpp fourtax.h partD.h Dfoil.h locusfile.h popZParent.h popZDstat.h popZpartD.h popZDfoil.h
	mpic++ -O3 -std=c++11  -Wall -c dstat_main.cpp

compDmpi : fourtax.o locusfile.o Stats.o DstatParent.o Dstat.o partD.o Dfoil.o popZParent.o popZDstat.o popZpartD.o popZDfoil.o dstat_main.o 
	mpic++ -O3 -std=c++11 -Wall -o compDmpi fourtax.o locusfile.o Stats.o DstatParent.o Dstat.o partD.o Dfoil.o popZParent.o popZDstat.o popZpartD.o popZDfoil.o dstat_main.o -lboost_program_options -lpthread

clean:
	rm *.o compDmpi

install:
	cp compDmpi /usr/local/bin/.
