#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "wallideal.h"
#include "lp.h"
#include "polyhedralcone.h"
#include "gfanapplication.h"
#include "polyhedralfan.h"
#include "halfopencone.h"
#include "matrix.h"
#include "regularsubdivision.h"
#include "triangulation2.h"
#include "log.h"

class RegularSubdivisionApplication : public GFanApplication
{
  SimpleOption conventionOption;
  SimpleOption triangulationOption;
  SimpleOption coordinateOption;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program takes a point configuration and a lifting vector a computes the corresponding regular subdivision PROJECTIVELY?.\n";
  }
  RegularSubdivisionApplication():
    conventionOption("--switch","Switch max min convention"),
    triangulationOption("--triangulation","Ensure that the output complex is a triangulation. This is done by perturbing the input vector lexicographically."),
    coordinateOption("--coordinate","Compute coordinatesof...")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_regularsubdivision";
  }

  int main()
  {
    FileParser P(Stdin);

    IntegerMatrix m=rowsToIntegerMatrix(P.parseIntegerVectorList());

    if(coordinateOption.getValue())
      {
      IntegerVectorList out;
      IntegerVectorList w=P.parseIntegerVectorList();
      m=m.transposed();
      Triangulation2 t(m);
      t.triangulate();
      for(IntegerVectorList::const_iterator i=w.begin();i!=w.end();i++)
        {
          WeightTermOrder T(*i);
          t.changeToTriangulationInducedBy(T);
          out.push_back(t.DFSResultantCoordinate());
        }
      pout<<out;
      return 0;
      }
    IntegerVector w=P.parseIntegerVector();

    if(conventionOption.getValue())w=-w;

    if(triangulationOption.getValue())
      {
        m=m.transposed();
        Triangulation2 t(m);
        t.triangulate();
        WeightTermOrder T(w);
        t.changeToTriangulationInducedBy(T);
          t.print(pout);
      }
    else
      printSetSetInt(stdout,regularSubdivision(m,w));

    return 0;
  }
};

static RegularSubdivisionApplication theApplication;
