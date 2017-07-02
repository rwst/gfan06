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
#include "symmetry.h"
#include "halfopencone.h"
#include "log.h"
#include "field_rationals.h"

class TropicalFunctionApplication : public GFanApplication
{
  SimpleOption exponentOption;
public:
  const char *helpText()
  {
    return "This program takes a polynomial and tropicalizes it. The output is piecewise linear function represented by a fan whose cones are the linear regions. Each ray of the fan gets the value of the tropical function assigned to it. In other words this program computes the normal fan of the Newton polytope of the input polynomial with additional information.";
  }
  TropicalFunctionApplication():
    exponentOption("--exponents","Tell program to read a list of exponent vectors instead.")
  {
    registerOptions();
  }
  const char *name()
  {
    return "_tropicalfunction";
  }
  void inner(PolynomialSet const &f)
  {
    PolyhedralFan F=PolyhedralFan::normalFanOfNewtonPolytope(*f.begin());

    {
      AsciiPrinter p(Stdout);
      PolyhedralFan a=F;
      a.printWithIndices(&p,FPF_default|FPF_values);
    }
  }
  int main()
  {
    FileParser P(Stdin);


    if(!exponentOption.getValue())
      {
        PolynomialSet f=P.parsePolynomialSetWithRing();
        inner(f);
      }
    else
      {
        IntegerVectorList exponents=P.parseIntegerVectorList();
        assert(exponents.size());
        int n=exponents.begin()->size();
        PolynomialRing R(Q,n);
        Polynomial p(R);
        for(IntegerVectorList::const_iterator i=exponents.begin();i!=exponents.end();i++)
          {
            p+=Term(Q.zHomomorphism(1),Monomial(R,*i));
          }
        PolynomialSet f(R);
        f.push_back(p);
        inner(f);
      }

    return 0;
  }
};

static TropicalFunctionApplication theApplication;
