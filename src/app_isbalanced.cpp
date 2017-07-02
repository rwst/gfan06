/*
 * app_isbalanced.cpp
 *
 *  Created on: May 25, 2011
 *      Author: anders
 */

#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "lp.h"
#include "gfanapplication.h"
#include "polyhedralcone.h"
#include "polyhedralfan.h"
#include "symmetry.h"
#include "log.h"
#include "matrix.h"
#include "linalg.h"
#include "polymakefile.h"

class BalancedApplication : public GFanApplication
{
  StringOption inputOption;
//  SimpleOption optionSymmetry;
  SimpleOption makeBalancedOption;
  SimpleOption facetComplexOption;
  SimpleOption coneOption;
public:
  const char *helpText()
  {
    return "This program checks if a fan is balanced\n";
  }
  BalancedApplication():
    inputOption("-i","Specify the name of the input file.","polymake.out"),
    makeBalancedOption("--makeBalanced","Assign positive multiplicities to cones to make the fan balanced instead. Fan is assumed to be pure."),
//    optionSymmetry("--symmetry","Reads in a fan stored with symmetry. The generators of the symmetry group must be given on the standard input.\n")
  facetComplexOption("--facetComplex","When making the fan balanced, first take facet complex"),
  coneOption("--cone","Output cone of all balanced weight functions.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_fanisbalanced";
  }

  int main()
  {
	    if(coneOption.getValue())
	      {
	        PolyhedralFan f=PolyhedralFan::readFan(inputOption.getValue(),true,0,0,0);
	        if(facetComplexOption.getValue())f=f.facetComplex();
	        int n=f.getAmbientDimension();

	        IntegerMatrix equations=f.balancingEquations();
	        int d=f.getMaxDimension();
	        int numberOfFacets=0;
	        for(PolyhedralFan::coneIterator i=f.conesBegin();i!=f.conesEnd();i++)numberOfFacets+=(i->dimension()==d);

	        PolyhedralCone multiplicityCone(IntegerVectorList(),equations.getRows(),numberOfFacets/*+ridges.size()*(d-1)*/);

//	        multiplicityCone=intersection(multiplicityCone,/*product(*/PolyhedralCone::positiveOrthant(numberOfFacets)/*,PolyhedralCone(ridges.size()*(d-1)))*/);

	        multiplicityCone.printAsFan(&pout);
	        return 0;
	      }



	        	if(makeBalancedOption.getValue())
      {
        PolyhedralFan f=PolyhedralFan::readFan(inputOption.getValue(),true,0,0,0);
        if(facetComplexOption.getValue())f=f.facetComplex();
        int n=f.getAmbientDimension();

        IntegerMatrix equations=f.balancingEquations();


#if 0
        {
          IntegerMatrix temp=equations.transposed();
          temp.prependRow(IntegerVector/*::allOnes*/(equations.getHeight()));
          IntegerMatrix temp2=IntegerMatrix::identity(numberOfFacets);
          temp2.prependRow(IntegerVector::allOnes(numberOfFacets));
          bool temp3=hasHomogeneousSolution(numberOfFacets+1,temp2.transposed().getRows(),temp.transposed().getRows());
          debug<<"hasHomogeneousSolution"<<temp3<<"\n";
          if(temp3)
            {
              debug<<"SKIPPING\n";
              return 0;
            }
        }
#endif

        int d=f.getMaxDimension();
        int numberOfFacets=0;
        for(PolyhedralFan::coneIterator i=f.conesBegin();i!=f.conesEnd();i++)numberOfFacets+=(i->dimension()==d);

        PolyhedralCone multiplicityCone(IntegerVectorList(),equations.getRows(),numberOfFacets/*+ridges.size()*(d-1)*/);

        multiplicityCone=intersection(multiplicityCone,/*product(*/PolyhedralCone::positiveOrthant(numberOfFacets)/*,PolyhedralCone(ridges.size()*(d-1)))*/);

        IntegerVector solution=multiplicityCone.getRelativeInteriorPoint().subvector(0,numberOfFacets);

        debug<<"SOLUTION"<<solution<<"\n";
        if(solution.supportAsZeroOneVector().sum()!=solution.size())
          {
            pout<<"Cannot be balanced\n";
            debug<<"Cannot be balanced\n";
            return 0;
          }
        solution=solution/gcdOfVector(solution);
        PolyhedralFan f2(n);
        debug<<"SOLUTION"<<solution<<"\n";
        int I=0;
        for(PolyhedralFan::coneIterator i=f.conesBegin();i!=f.conesEnd();i++,I++)
          {
            PolyhedralCone temp=*i;
            temp.setMultiplicity(solution[I]);
            f2.insert(temp);
          }
        SymmetryGroup S(n);
        f2.printWithIndices(&pout,
                            FPF_multiplicities|FPF_default |
                            /*(optionSymmetry.getValue()?FPF_group|FPF_conesCompressed:0)|*/
                            /*(optionIgnoreCones.getValue()?0:FPF_conesExpanded)|*/
                            FPF_maximalCones|FPF_cones,
                            &S);

//        pout<<multiplicityCone.getRelativeInteriorPoint();
        return 0;
      }



    /*    if(optionSymmetry.getValue())
      {
        IntegerVectorList generators=FileParser(Stdin).parseIntegerVectorList();
        s.computeClosure(generators);
      }*/
    PolyhedralFan f=PolyhedralFan::readFan(inputOption.getValue(),true,0,0,0);

    PolyhedralFan skeleton=f.facetComplex();
    bool isBalanced=true;
    for(PolyhedralFan::coneIterator i=skeleton.conesBegin();i!=skeleton.conesEnd();i++)
      {
        log1 debug<<"checking"
            "\n";
        PolyhedralFan L=f.link(i->getRelativeInteriorPoint());
        int d=L.getMaxDimension();
        IntegerVector sum(L.getAmbientDimension());
        for(PolyhedralFan::coneIterator j=L.conesBegin();j!=L.conesEnd();j++)
          sum=sum+(j->getMultiplicity()*(j->semiGroupGeneratorOfRay()));
        if(!L.conesBegin()->linealitySpace().contains(sum))
          {
            log1 debug<<"Not balanced at:"<<i->getRelativeInteriorPoint()<<"\n";
            log1 for(PolyhedralFan::coneIterator j=L.conesBegin();j!=L.conesEnd();j++)
              debug<<j->getMultiplicity()<<(j->semiGroupGeneratorOfRay())<<"\n";
            log1 debug<<"Sum"<<sum;
            isBalanced=false;
            break;
          }
      }
    pout<<isBalanced<<"\n";
    return 0;
  }
};

static BalancedApplication theApplication;

