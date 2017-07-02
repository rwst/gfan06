#include <iostream>
#include <stdlib.h>
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
#include "polyhedralcone.h"
#include "gfanapplication.h"
#include "saturation.h"
#include "field_rationals.h"
#include "field_zmodpz.h"
#include "field_rationalfunctions.h"
#include "symmetry.h"
#include "linalg.h"
#include "fieldlp.h"
#include "integer.h"
#include "polynomialgcd.h"
#include "packedmonomial.h"
#include "gfanlib_zcone.h"

using namespace gfan;

class TropicalVarietySpanApplication : public GFanApplication
{
public:
	bool includeInDefaultInstallation() // Not included since the program does not follow the specifications of the help text. Moreover, doing this using Hadamard products is extremely slow. This should be replaced by the algorithm in the paper by Kahle, Katthan and Jensen.
	{
		return false;
	}
  const char *helpText()
  {
    return "This program takes a ring and an ideal as input and computes the span of the tropical variety defined by the ideal.\n";
  }
  TropicalVarietySpanApplication()
  {
    registerOptions();
  }

  const char *name()
  {
    return "_tropicalvarietyspan";
  }

  PolynomialSet hadamardProduct(PolynomialSet const &I, PolynomialSet const &J)
  {
	  PolynomialRing R=I.getRing();
	  int n=R.getNumberOfVariables();
	  PolynomialRing bigRing(R.getField(),n*3);
	  vector<int> var1;
	  vector<int> var2;
	  for(int i=0;i<n;i++)
	  {
		  var1.push_back(i+n);
		  var2.push_back(i+2*n);
	  }
	  PolynomialSet all(bigRing);
	  all.unionSet(I.embeddedInto2(bigRing,var1));
	  all.unionSet(J.embeddedInto2(bigRing,var2));
	  for(int i=0;i<n;i++)
		  all.push_back(bigRing.ithVariable(i)-bigRing.ithVariable(n+i)*bigRing.ithVariable(2*n+i));

	  IntegerVectorList l;
//	  l.push_back();
	  l.push_back(IntegerVector(3*n));
	  l.push_back(concatenation(IntegerVector(n+n),IntegerVector::allOnes(n)));
	  MatrixTermOrder T(l);
	  buchberger(&all,T,true);
	  return all.polynomialRingIntersection(R);
  }

  int main()
  {
	  FileParser P(Stdin);

	  PolynomialSet g=P.parsePolynomialSetWithRing();

	  PolynomialSet s=g;
	  s=hadamardProduct(s,g);
	  pout<<s;
	  s=hadamardProduct(s,s);
	  pout<<s;


//	  pout<<hadamardProduct(g,g);

	  return 0;
  }
};

static TropicalVarietySpanApplication theApplication;
