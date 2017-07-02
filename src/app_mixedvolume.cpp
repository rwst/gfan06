#include <iostream>
#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "matrix.h"
#include "gfanlib_matrix.h"
//#include "mixedvolume.h"
#include "gfanapplication.h"
#include "gfanlib_mixedvolume.h"
#include "log.h"

using namespace gfan;
using namespace gfan::MixedVolumeExamples;

class MixedVolumeApplication : public GFanApplication
{
	SimpleOption optionVectorInput;
	IntegerOption optionCyclic;
	IntegerOption optionNoon;
	IntegerOption optionChandra;
	IntegerOption optionKatsura;
	IntegerOption optionGaukwa;
	IntegerOption optionEco;
	IntegerOption optionNThreads;
	IntegerOption optionSteps;
public:
  const char *helpText()
  {
    return "This program computes the mixed volume of the Newton polytopes of a list of polynomials. The ring is specified on the input. After this follows the list of polynomials.\n";
  }
  MixedVolumeApplication():
	  optionVectorInput("--vectorinput","Read in a list of point configurations instead of a polynomial ring and a list of polynomials."),
	  optionCyclic("--cyclic","Use cyclic-n example instead of reading input."),
	  optionNoon("--noon","Use Noonburg-n example instead of reading input."),
	  optionChandra("--chandra","Use Chandrasekhar-n example instead of reading input."),
	  optionKatsura("--katsura","Use Katsura-n example instead of reading input."),/* Note that Verschelde's mixed volumes for the Katsura examples do not match those produced by gfan. The configurations do not seem to be identical for all n. */
	  optionGaukwa("--gaukwa","Use Gaukwa-n example instead of reading input."),
	  optionEco("--eco","Use Eco-n example instead of reading input."),
	  optionNThreads("-j","Number of threads"),
	  optionSteps("-s","Number of steps", 500)
  {
	  optionSteps.hide();
	  registerOptions();
  }
  const char *name()
  {
    return "_mixedvolume";
  }
  IntMatrix rowsToIntegerMatrix(IntegerVectorList const &l)
  {
	  assert(l.size());
	  IntMatrix ret(l.size(),l.front().size());
	  IntegerVectorList::const_iterator I=l.begin();
	  for(int i=0;i<ret.getHeight();i++,I++)
		  for(int j=0;j<ret.getWidth();j++)
			  ret[i][j]=(*I)[j];
	  return ret;
  }
  int main()
  {
		int nthreads=optionNThreads.getValue();
		int steps=optionSteps.getValue();
		vector<IntMatrix> tuple;

		if(optionCyclic.getValue())
			tuple=cyclic(optionCyclic.getValue());
		else if(optionNoon.getValue())
			tuple=noon(optionNoon.getValue());
		else if(optionChandra.getValue())
			tuple=chandra(optionChandra.getValue());
		else if(optionKatsura.getValue())
			tuple=katsura(optionKatsura.getValue());
		else if(optionGaukwa.getValue())
			tuple=gaukwa(optionGaukwa.getValue());
		else if(optionEco.getValue())
			tuple=eco(optionEco.getValue());
		else
		{
			IntegerVectorListList tuple1;
			if(optionVectorInput.getValue())
			{
				tuple1=FileParser(Stdin).parseIntegerVectorListList();
			}
			else
			{
				PolynomialSet g=FileParser(Stdin).parsePolynomialSetWithRing();

				tuple1=g.exponents();
			}
			for(IntegerVectorListList::const_iterator i=tuple1.begin();i!=tuple1.end();i++)
				tuple.push_back(rowsToIntegerMatrix(*i).transposed());
		}

		log1 for(auto i=tuple.begin();i!=tuple.end();i++)std::cerr<<*i;

		cout<<mixedVolume(tuple,nthreads,steps)<<"\n";
		//    cout << mixedVolume(s) << endl;

		return 0;
  }
};

static MixedVolumeApplication theApplication;
