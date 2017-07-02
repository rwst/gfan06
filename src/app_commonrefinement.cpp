#include "parser.h"
#include "printer.h"
#include "lp.h"
#include "gfanapplication.h"
#include "polyhedralcone.h"
#include "polyhedralfan.h"
#include "polymakefile.h"

class CommonRefinementApplication : public GFanApplication
{
  StringOption input1Option;
  StringOption input2Option;
  SimpleOption stableOption;
public:
  const char *helpText()
  {
    return "This program takes two polyhedral fans and computes their common refinement.\n";
  }
  CommonRefinementApplication():
    input1Option("-i1","Specify the name of the first input file.","polymake.out"),
    input2Option("-i2","Specify the name of the second input file.","polymake.out"),
    stableOption("--stable","Compute the stable intersection.")
  {
    //    stableOption.hide();
    registerOptions();
  }

  const char *name()
  {
    return "_fancommonrefinement";
  }

  int main()
  {
    PolyhedralFan f1=PolyhedralFan::readFan(input1Option.getValue());
    PolyhedralFan f2=PolyhedralFan::readFan(input2Option.getValue());

    PolyhedralFan f=refinement(f1,f2,-1,false,stableOption.getValue());

    AsciiPrinter P(Stdout);

    f.printWithIndices(&P,FPF_default/*|(stableOption.getValue()?FPF_multiplicities:0)|FPF_values*/);
    /* TODO: If one wants to do intersection theory as Fulton and Sturmfels
     * it would be very convenient that multiplicities are output here.
     * Unfortunately the refinement command does not compute multiplicities at the moment.
     * This should be fixed to make gfan more useful.
     * */


    return 0;
  }
};

static CommonRefinementApplication theApplication;
