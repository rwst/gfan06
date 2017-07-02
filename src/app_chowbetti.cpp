/*
 *  app_chowbetti.cpp
 *
 *  Created on: June 21, 2016
 *
 */

#include "gfanapplication.h"
#include "polyhedralcone.h"
#include "polyhedralfan.h"
#include "log.h"
#include "matrix.h"
#include "linalg.h"
#include "polymakefile.h"

class ChowBettiApplication : public GFanApplication
{
  StringOption inputOption;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program computes the Chow Betti numbers of a pure fan, i.e. dimensions of the spaces of tropical cycles of each dimension supported on the input fan. \n";
  }
  ChowBettiApplication():
    inputOption("-i","Specify the name of the input file.","polymake.out")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_chowbetti";
  }

  int main()
  {
    PolyhedralFan f=PolyhedralFan::readFan(inputOption.getValue(),true,0,0,0);
    int n=f.getMaxDimension();
    IntegerVector ChowBetti(n);
    // Now we are assuming that the input fan is pure.
    // The Chow Betti numbers make sense for non-pure fans as well.

    while(f.getMaxDimension() != 0)
    {
      pout << "dim = " << f.getMaxDimension() << "\n";
      FieldMatrix equations = integerMatrixToFieldMatrix(f.balancingEquations(),Q);
      ChowBetti[n-f.getMaxDimension()] = equations.getWidth() - equations.reduceAndComputeRank();
      f=f.facetComplex();
      pout <<  ChowBetti << "\n";
    }
    pout <<  ChowBetti << "\n";
    return 0;
  }
};

static ChowBettiApplication theApplication;
