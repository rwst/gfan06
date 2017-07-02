/*
 * app_tropicalcurve.cpp
 *
 *  Created on: Apr 23, 2014
 *      Author: anders
 */

#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "reversesearch.h"
#include "termorder.h"
#include "ep_standard.h"
#include "ep_xfig.h"
#include "gfanapplication.h"
#include "polyhedralcone.h"
#include "polyhedralfan.h"
#include "tropical.h"
#include "tropical2.h"
#include "tropicalbasis.h"
#include "tropicalcurve.h"
#include "dimension.h"
#include "field_rationalfunctions2.h"
#include "log.h"

class TropicalCurveApplication : public GFanApplication
{

//  SimpleOption optionHomogenize;
	SimpleOption optionSingleRay;
	IntegerOption optionParameters;
public:
  const char *helpText()
  {
    return "This program computes a tropical basis for an ideal defining a tropical curve. Defining a tropical curve means that the Krull dimension of R/I is at most 1 + the dimension of the homogeneity space of I where R is the polynomial ring. The input is a generating set for the ideal. If the input is not homogeneous option -h must be used.\n";
  }
  TropicalCurveApplication():
//    optionHomogenize("-h","Homogenise the input before computing a tropical basis and dehomogenise the output. This is needed if the input generators are not already homogeneous.")
	  optionParameters("--parameters","With this option you can specify how many variables to treat as parameters instead of variables. This makes it possible to do computations where the coefficient field is the field of rational functions in the parameters.",0),
	  optionSingleRay("--singleray","Only compute a single ray of the curve.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_tropicalcurve";
  }

  IntegerVectorList tropicalCurveWithPreprocessing(PolynomialSet const &theInput)
  {
	  {
		  list<int> temp=theInput.multiDeHomogenizationToKeep();
		  for(list<int>::const_iterator i=temp.begin();i!=temp.end();i++)debug<<*i<<"\n";
	  }
	  debug<<theInput.getRing()<<theInput;
	  debug<<"A\n";
	  PolynomialRing R=theInput.getRing();
	  PolynomialSet groebnerBasis=theInput;
	  debug<<"B-----\n"<<R<<"\n"<<groebnerBasis;
	  buchberger(&groebnerBasis,WeightReverseLexicographicTermOrder(IntegerVector::allOnes(R.getNumberOfVariables())),true);


	  debug<<"B-----\n"<<R<<"\n"<<groebnerBasis;
//	  assert(0);
debug<<"C-----\n";

	  list<int> toKeep=groebnerBasis.multiDeHomogenizationToKeep();
	  PolynomialSet g1=groebnerBasis.multiDeHomogenization();
		PolynomialSet g=g1.homogenization(g1.getRing().withVariablesAppended("H"));
	//	saturatedIdeal(g);
		  buchberger(&g,WeightReverseLexicographicTermOrder(IntegerVector::allOnes(g.getRing().getNumberOfVariables())));

		int d=krullDimension(g);
		debug<<"d="<<d<<"\n";
		assert(d==1||d==2);
		if(d==1)return IntegerVectorList();

		debug<<g;
		PolynomialSet dehom=g.deHomogenization();

		debug<<dehom;
		IntegerVectorList rays=tropicalCurve(dehom,optionSingleRay.getValue());

		debug<< "Reduced coordinates"<<rays<<"\n";
/*		IntegerVectorList rays2;
		for(IntegerVectorList::const_iterator i=rays.begin();i!=rays.end();i++)
		{
			rays2.push_back(i->subvector(0,i->size()-1)-(*i)[i->size()-1]*IntegerVector::allOnes(i->size()-1));
		}
		debug<< "Dehomogenized"<<rays2<<"\n";
*/

		debug<<toKeep;
		IntegerVectorList rays2;
		for(IntegerVectorList::const_iterator i=rays.begin();i!=rays.end();i++)
		{
			IntegerVector v(theInput.getRing().getNumberOfVariables());
//			debug<<v;
			int J=0;
			for(list<int>::const_iterator j=toKeep.begin();j!=toKeep.end();j++,J++)
				v[*j]=(*i)[J];
			rays2.push_back(v);
		}
	    return rays2;
  }
  int main()
  {
    FileParser P(Stdin);

    PolynomialSet theInput=P.parsePolynomialSetWithRing();
    if(optionParameters.getValue())
    	theInput=makeVariablesParameters(
    			makeVariablesParameters(theInput.getRing(),optionParameters.getValue())
    			,theInput);

    AsciiPrinter(Stdout).printVectorList(tropicalCurveWithPreprocessing(theInput));

    return 0;
  }
};

static TropicalCurveApplication theApplication;
