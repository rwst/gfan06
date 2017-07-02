/*
 * app_polytopealgebra.cpp
 *
 *  Created on: Nov 30, 2011
 *      Author: anders
 */

#include <iostream>
#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "reversesearch.h"
#include "breadthfirstsearch.h"
#include "termorder.h"
#include "ep_standard.h"
#include "ep_xfig.h"
#include "gfanapplication.h"
#include "polyhedralfan.h"
#include "minkowskidual.h"
#include "log.h"
#include "tropical_weildivisor.h"
#include "field_rationals.h"
#include "linalg.h"

class PolytopeAlgebraApplication : public GFanApplication
{
public:
  bool includeInDefaultInstallation() // Not included since the program has not been documented and is likely not working.
  {
	return false;
  }
  const char *helpText()
  {
    return "This does some experimental computations in the polytope algebra.\n";
  }
  PolytopeAlgebraApplication()
  {
    registerOptions();
  }

  const char *name()
  {
    return "_polytopealgebra";
  }

  PolyhedralFan simpleRationalFunction(PolyhedralFan const &support, IntegerVector const &withValueOne)
  {
//D(withValueOne);
    int n=support.getAmbientDimension();
    PolyhedralFan ret(n);
    for(PolyhedralFan::coneIterator i=support.conesBegin();i!=support.conesEnd();i++)
      if(!i->contains(withValueOne))
        {
          PolyhedralCone temp(*i);
          temp.setLinearForm(IntegerVector(n));
          ret.insert(temp);
        }
      else
        {
          PolyhedralCone temp(*i);
          IntegerVectorList rays=i->extremeRays();
          IntegerVectorList lines=i->generatorsOfLinealitySpace();

          FieldMatrix LFA(Q,lines.size()+rays.size(),n);
          FieldVector LFB(Q,lines.size()+rays.size()+n);
          int J=0;
          for(IntegerVectorList::const_iterator j=rays.begin();j!=rays.end();j++,J++)
              {
              LFA[J]=integerVectorToFieldVector(*j,Q);
              LFB[J]=Q.zHomomorphism((*j-withValueOne).isZero()*120);
              }
    //      D(LFA);
    //      D(LFB);
          FieldVector LFX=LFA.solver().canonicalize(LFB);
          if(LFX.subvector(0,LFX.size()-n).isZero())
            {
              temp.setLinearForm(fieldVectorToIntegerVector(LFX.subvector(LFX.size()-n,LFX.size())));
      //        D(temp.getLinearForm());
            }
          else
            {
              cerr<<"Values on cone are not linear" <<endl;
              assert(0);
            }


          ret.insert(temp);
        }
    return ret;
  }

  int main()
  {
    LexicographicTermOrder myOrder;

    PolynomialSet g=FileParser(Stdin).parsePolynomialSetWithRing();
    g.sort_();
    g.markAndScale(myOrder);

    int k=2;

    PolyhedralFan supportFan=PolyhedralFan::normalFanOfNewtonPolytope(*g.begin());

    supportFan=supportFan.triangulation();
    supportFan.printWithIndices(&debug,FPF_multiplicities|FPF_default);

    PolyhedralFan a=supportFan;

    for(int i=0;i<k;i++)
      a=a.facetComplex();

    IntegerMatrix E=a.balancingEquations();
    pout<<E.getRows();
    pout<<"HEIGHT "<<E.getHeight()<<"WIDTH "<<E.getWidth()<<"RANK "<<::rank(E)<<"\n";

    IntegerVectorList rays=a.getRaysInPrintingOrder(0);
    pout<<"RAYS"<<rays;

    IntegerVectorList generators;
    int n=supportFan.getAmbientDimension();
    for(IntegerVectorList::const_iterator i=rays.begin();i!=rays.end();i++)
      for(IntegerVectorList::const_iterator j=rays.begin();j!=rays.end();j++)
        {
          PolyhedralFan A=simpleRationalFunction(supportFan,*i);
         // A.printWithIndices(&debug,FPF_multiplicities|FPF_default|FPF_values);
          PolyhedralFan B=simpleRationalFunction(supportFan,*j);

          PolyhedralFan product=PolyhedralFan::fullSpace(n);
          product=weilDivisor(product,A);
          product=weilDivisor(product,B);


    //      debug>>
    //      product.printWithIndices(&debug,FPF_multiplicities|FPF_default);
          IntegerVector generator(E.getWidth());
          int K=0;
          for(PolyhedralFan::coneIterator k=a.conesBegin();k!=a.conesEnd();k++,K++)
            {
              if(product.containsInSupport(k->getRelativeInteriorPoint()))
                generator[K]=product.coneContaining(k->getRelativeInteriorPoint()).getMultiplicity();
            }
          generators.push_back(generator);
        }
    IntegerMatrix F=rowsToIntegerMatrix(generators,E.getWidth());
    pout<<F.getRows();
    pout<<"HEIGHT "<<F.getHeight()<<"WIDTH "<<F.getWidth()<<"RANK "<<::rank(F)<<"\n";

    pout<<"PRODUCT "<<(E*(F.transposed())).getRows();

    return 0;
  }
};

static PolytopeAlgebraApplication theApplication;

