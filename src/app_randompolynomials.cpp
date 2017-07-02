/*
 * app_randompolynomials.cpp
 *
 *  Created on: Feb 2, 2012
 *      Author: anders
 */

#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "minors.h"
#include "field_rationals.h"
#include <iostream>

class RandomPolynomialsApplication : public GFanApplication
{
  IntegerOption nvarsOption;
  IntegerOption npolysOption;
  IntegerOption ntermsOption;
  IntegerOption maxDegreeOption;
public:
  bool includeInDefaultInstallation() // Not included since we don't want to support this program in the future.
  {
    return false;
  }
  const char *helpText()
  {
    return "This program produces a ring with random polynomials.\n";
  }
  RandomPolynomialsApplication():
    nvarsOption("--nvars","Specify number of variables.",1),
    npolysOption("--npolys","Specify number of polynomials.",1),
    ntermsOption("--nterms","Specify number of terms per poly.",1),
    maxDegreeOption("--maxdegree","Set maximal degree in each variable.",1)
  {
    registerOptions();
  }
  const char *name()
  {
    return "_randompolynomials";
  }

  int main()
  {

    PolynomialRing R(Q,nvarsOption.getValue());

    PolynomialSet g(R);

    srand(time(NULL));

    for(int i=0;i<npolysOption.getValue();i++)
    {
    	Polynomial p(R);

    	for(int j=0;j<ntermsOption.getValue();j++)
    	{
    		IntegerVector e(nvarsOption.getValue());
//    		cerr<<maxDegreeOption.getValue();
    		for(int k=0;k<e.size();k++)e[k]=(rand())%((unsigned int)maxDegreeOption.getValue());
    		p+=Term(Q.zHomomorphism(1),Monomial(R,e));
    	}

    	g.push_back(p);
    }

//		  p+=Term(R2.getField().zHomomorphism(1),Monomial(R2,e));

    pout << R<<g;

    return 0;
  }
};

static RandomPolynomialsApplication theApplication;
