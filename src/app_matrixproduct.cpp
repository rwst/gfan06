/*
 * app_matrixproduct.cpp
 *
 *  Created on: Jan 4, 2011
 *      Author: anders
 */

#include "parser.h"
#include "printer.h"
#include "gfanapplication.h"
#include "matrix.h"

class MatrixProductApplication : public GFanApplication
{
public:
  SimpleOption optionTropical;
  SimpleOption optionNegate;
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program computes the product of two matrices.\n";
  }
  MatrixProductApplication():
    optionTropical("--tropical","Do the computation in the max-plus semi-ring."),
    optionNegate("--negate","Change the sign of a matrix instead.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_matrixproduct";
  }

  int main()
  {
    FileParser P(Stdin);

    IntegerMatrix A=rowsToIntegerMatrix(P.parseIntegerVectorList());
    if(optionNegate.getValue())
      pout<<(-A).getRows();
    else
      {
        IntegerMatrix B=rowsToIntegerMatrix(P.parseIntegerVectorList());

        pout<<((optionTropical.getValue())?tropicalProduct(A,B):A*B).getRows();
      }

    return 0;
  }
};

static MatrixProductApplication theApplication;
