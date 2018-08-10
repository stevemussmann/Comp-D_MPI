AM_CPPFLAGS = -Wall -std=c++11 -O3 $(BOOST_CPPFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS)

#this lists the binaries to produce
bin_PROGRAMS = compDmpi fn
compDmpi_SOURCES = dstat_main.cpp fourtax.h partD.h Dfoil.h locusfile.h popZParent.h popZDstat.h popZpartD.h popZDfoil.h
fn_SOURCES = fn_main.cpp

noinst_LIBRARIES = libstats.a libdstatparent.a liblocusfile.a libfourtax.a libdstat.a libpartd.a libdfoil.a libpopzparent.a libpopzdstat.a libpopzpartd.a libpopzdfoil.a libfnfiles.a libfnstats.a

#compd
libstats_a_SOURCES = Stats.h Stats.cpp

libdstatparent_a_SOURCES = DstatParent.h DstatParent.cpp Stats.h

liblocusfile_a_SOURCES = locusfile.h locusfile.cpp

libfourtax_a_SOURCES = fourtax.h fourtax.cpp

libdstat_a_SOURCES = Dstat.h Dstat.cpp fourtax.h DstatParent.h Stats.h

libpartd_a_SOURCES = partD.h partD.cpp fourtax.h DstatParent.h Stats.h

libdfoil_a_SOURCES = Dfoil.h Dfoil.cpp fourtax.h DstatParent.h Stats.h

libpopzparent_a_SOURCES = popZParent.h popZParent.cpp Stats.h Dstat.h partD.h Dfoil.h

libpopzdstat_a_SOURCES = popZDstat.h popZDstat.cpp Stats.h Dstat.h DstatParent.h

libpopzpartd_a_SOURCES = popZpartD.h popZpartD.cpp Stats.h Dstat.h DstatParent.h

libpopzdfoil_a_SOURCES = popZDfoil.h popZDfoil.cpp Stats.h Dstat.h DstatParent.h

compDmpi_LDFLAGS = $(MPI_CXXLDFLAGS) -lboost_program_options -lpthread
compDmpi_CPPFLAGS = -Wall -std=c++11 -O3
compDmpi_LDADD = fourtax.o locusfile.o Stats.o DstatParent.o Dstat.o partD.o Dfoil.o popZParent.o popZDstat.o popZpartD.o popZDfoil.o -lboost_program_options -lpthread

#fn
libfnfiles_a_SOURCES = fnFiles.cpp fnFiles.h

libfnstats_a_SOURCES = fnStats.cpp fnStats.h

fn_LDFLAGS = 
fn_CPPFLAGS = -Wall -std=c++11 -O3
fn_LDADD = libfnfiles.a libfnstats.a -lboost_program_options

