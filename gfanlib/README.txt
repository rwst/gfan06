This is gfanlib for doing polyhedral computations related to polynomial rings.

Contributors are:
	     Anders Nedergaard Jensen (main contributor)
	     Bjarne Knudsen (parallel acyclic graph traverser)
	     Yue Ren (adjustments of the library to fit Singular)
	     The Singular Team (Makefile.in and configure.ac)

Gfanlib relies on cddlib by Komei Fukuda and on gmp.

Gfanlib has two major feature:

1) high-level exact polyhedral cone and polyhedral fan classes
2) fast exact mixed volume computation for lattice polytopes with overflow checking

In particular, gfanlib is missing the Groebner basis part of gfan.

If only mixed volume computations are desired, then cddlib is not used and gmp is not essential and these libraries can be removed with a little work.

To compile on linux-like systems do:

export CPPFLAGS="-I/usr/include/cdd" 
./configure
make

To obtain speed for mixed volume computations as advertised elsewhere, change the optimisation flags 
./configure CXXFLAGS='-O3 -mavx -msse2 -finline-limit=1000'
before making and use gcc version >=4.7.2.

To demonstrate usage of the library, we provide the following session with
the CLING C++ interpreter:

anders@isfrun ~/gfan-svn/jensen-gfan/trunk/gfanlib $ ~/software/cling3/cling-Ubuntu-12.04-64bit-fca52f08e2/bin/cling -l /usr/lib/libgmp.so.3.5.2  -DGMPRATIONAL -std=c++0x

****************** CLING ******************
* Type C++ code and press enter to run it *
*             Type .q to exit             *
*******************************************
[cling]$ #include "gfanlib_circuittableint.cpp"
[cling]$ #include "gfanlib_paralleltraverser.cpp"
[cling]$ #include "gfanlib_mixedvolume.cpp"
[cling]$ using namespace std;
[cling]$ using namespace gfan;
[cling]$ using namespace gfan::MixedVolumeExamples;
[cling]$ auto s={cyclic(2),cyclic(3),cyclic(5),cyclic(7)};
[cling]$ for(auto v:s)cout<<mixedVolume(v)<<" ";cout<<endl;
2 6 70 924 
[cling]$ .q
anders@isfrun ~/gfan-svn/jensen-gfan/trunk/gfanlib $ 

And a compiled example with 1 and 16 threads:

#include <iostream>
#include "gfanlib.h"
using namespace gfan;
using namespace std;
int main()
{
    try{
	cout<<mixedVolume(MixedVolumeExamples::cyclic(12))<<endl;
	cout<<mixedVolume(MixedVolumeExamples::cyclic(12),16)<<endl;
	cout<<mixedVolume(MixedVolumeExamples::gaukwa(7),16)<<endl;
    }
    catch (...)
    {
	cerr<<"Error - most likely an integer overflow."<<endl;
    }
    return 0;
}

to be compiled with 
g++ -std=c++0x test.cpp -lgmp -lpthread libgfan.a 
