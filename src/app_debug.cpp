#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "groebnerengine.h"
#include "termorder.h"
#include "gfanapplication.h"
#include "log.h"

class DebugApplication : public GFanApplication
{
public:
	bool includeInDefaultInstallation()
	{
		return false;
	}
	const char *helpText()
	{
		return "This program is just for testing and debugging.\n";
	}
  DebugApplication()
  {
    registerOptions();
  }

  const char *name()
  {
    return "_debug";
  }

  int main()
  {
#if 1
	  PolynomialSet g=StringParser("{e^2*g^2-3*d*e^2*i,"
    		"d*e*f^6*g^2*h-1/2*d*e*f^7*g*i-3*d^2*e*f^6*h*i+9/2*d^2*e*f^7*j,"
    		"d*e^2*f^4*g*i-9*d^2*e^2*f^4*j,"
    		"d*e^3*g*i-9*d^2*e^3*j,"
    		"d^2*e*f^4*g^2*h^2-2/9*d*e*f^5*g^3*h+1/9*d*e*f^6*g^2*i+1/3*d^2*e*f^6*i^2-2*d^2*e*f^6*g*j-3*d^3*e*f^4*h^2*i+6*d^3*e*f^5*h*j,"
    		"d^2*e*f^7*g*h*i-2*d^2*e*f^8*i^2+6*d^2*e*f^8*g*j-9*d^3*e*f^7*h*j,"
    		"d^2*e*f^8*g^2*i^2-4*d^2*e*f^8*g^3*j-4*d^3*e*f^8*i^3+18*d^3*e*f^8*g*i*j-27*d^4*e*f^8*j^2,"
    		"d^2*e^2*f*g*h*i-d^2*e^2*f^2*i^2+3*d^2*e^2*f^2*g*j-9*d^3*e^2*f*h*j,"
    		"d^2*e^2*f^4*i^2-3*d^2*e^2*f^4*g*j,"
    		"d^2*e^3*i^2-3*d^2*e^3*g*j,"
    		"d^3*f^7*g^2*h^4-8/9*d^2*f^8*g^3*h^3+4/3*d^2*f^9*g^2*h^2*i-2/3*d^2*f^10*g*h*i^2+1/9*d^2*f^11*i^3+2*d^3*f^8*g*h^3*i-2*d^3*f^9*h^2*i^2-6*d^3*f^9*g*h^2*j+6*d^3*f^10*h*i*j-3*d^3*f^11*j^2-3*d^4*f^7*h^4*i+6*d^4*f^8*h^3*j,"
    		"d^3*e*f*g^2*h^3+1/27*d*e*f^3*g^4*h-1/54*d*e*f^4*g^3*i-5/9*d^2*e*f^2*g^3*h^2+1/3*d^2*e*f^3*g^2*h*i-1/6*d^2*e*f^4*g*i^2+1/2*d^2*e*f^4*g^2*j+d^3*e*f^2*g*h^2*i-4*d^3*e*f^3*g*h*j+1/2*d^3*e*f^4*i*j-3*d^4*e*f*h^3*i+6*d^4*e*f^2*h^2*j,"
    		"d^3*e*f^7*h*i^2-1/2*d^2*e*f^8*g*i^2+2*d^2*e*f^8*g^2*j-3*d^3*e*f^7*g*h*j-3/2*d^3*e*f^8*i*j,"
    		"d^3*e^2*f*h*i^2-1/3*d^2*e^2*f^2*g*i^2-3*d^3*e^2*f*g*h*j+3*d^3*e^2*f^2*i*j}").parsePolynomialSetWithRing();
    WeightReverseLexicographicTermOrder T(StringParser("(0,0,0,0,0,1,0,0,0,0)").parseIntegerVector());
#else
	  PolynomialSet g=StringParser("{a+b+c,b+c}").parsePolynomialSetWithRing();
  WeightReverseLexicographicTermOrder T(StringParser("(100,10,0)").parseIntegerVector());
#endif
    pout<<g;

    pout<<GE_groebnerBasis(g,T,true,false);
    return 0;
  }
};

static DebugApplication theApplication;

