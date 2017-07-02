/*
 * app_resultantfan.cpp
 *
 *  Created on: Mar 17, 2011
 *      Author: anders
 */
#include <assert.h>
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
#include "timer.h"
#include "log.h"
#include "matrix.h"
#include "lll.h"
#include "polyhedralfan.h"
#include "linalg.h"
#include "determinant.h"
#include "triangulation.h"
#include "intsinpolytope.h"
#include "graph.h"
#include "halfopencone.h"
#include "myassert.h"
#include "triangulation2.h"
#include "bsptree.h"

#include "traverser_secondaryfan.h"
#include "traverser_resultantfan.h"
#include "traverser_resultantfanspecialization.h"
#include "symmetrictraversal.h"
#include "traverser_bsptree.h"

#include <iostream>
#include <algorithm>

class ResultantFanApplication : public GFanApplication
{
 // StringOption inputOption;
  SimpleOption oldOption;
  SimpleOption optionCodimension;
  SimpleOption symmetryOption;
  SimpleOption optionIgnoreCones;
  SimpleOption optionSpecialization;
  SimpleOption optionVectorInput;
  SimpleOption optionProjection;
  SimpleOption optionPrettyPrint;
public:
  const char *helpText()
  {
    return "This program computes the resultant fan as defined in \"Computing Tropical Resultants\" by Jensen and Yu. The input is a polynomial ring followed by polynomials, whose coefficients are ignored. The output is the fan of coefficients such that the input system has a tropical solution.\n";
  }
  ResultantFanApplication():
   // inputOption("-i","Specify the name of the input file.",0),
    symmetryOption("--symmetry","Tells the program to read in generators for a group of symmetries (subgroup of $S_n$) after having read in the vector configuration. The program DOES NOT checks that the configuration stays fixed when permuting the variables with respect to elements in the group. The output is grouped according to the symmetry.\n"),
    optionIgnoreCones("--nocones","Tells the program not to output the CONES and MAXIMAL_CONES sections, but still output CONES_COMPRESSED and MAXIMAL_CONES_COMPRESSED if --symmetry is used.\n"),
    oldOption("--old","Use old stupid algorithm"),
    optionCodimension("--codimension","Compute only the codimension of the resultant fan and return.\n"),
    optionSpecialization("--special","Read in a zero-one vector from the standard input and specialize all variables with a one. That is, compute the stable intersection of the resultant fan with the subspace where the variables with a one in the vector are forced to zero. AT THE MOMENT ALSO A RELATIVE INTERIOR POINT OF A STARTING CONE IS READ.\n"),
    optionVectorInput("--vectorinput","Read in a list of point configurations instead of a polynomial ring and a list of polynomials.\n"),
    optionProjection("--projection","Use the projection method to compute the resultant fan. This works only if the resultant fan is a hypersurface. If this option is combined with --special, then the output fan lives in the subspace of the non-specialized coordinates.\n"),
    optionPrettyPrint("--pretty","Pretty print point configuration in TeX math and exit.\n")
  {
	  oldOption.hide();
	  optionPrettyPrint.hide();
    registerOptions();
  }

  const char *name()
  {
    return "_resultantfan";
  }


    vector<pair<int,int> > tupleToIntervals(IntegerVectorListList const &tuple)
    {
      vector<pair<int,int> > ret;
      int last=0;
      for(IntegerVectorListList::const_iterator i=tuple.begin();i!=tuple.end();i++)
        {
          ret.push_back(pair<int,int>(last,last+i->size()));
          last+=i->size();
        }
      return ret;
    }
    int findD(IntegerVectorListList const &tuple)
    {
      int d=0;
      for(IntegerVectorListList::const_iterator i=tuple.begin();i!=tuple.end();i++)
        for(IntegerVectorList::const_iterator j=i->begin();j!=i->end();j++)
          {d=j->size();goto leave;}
      leave:
      for(IntegerVectorListList::const_iterator i=tuple.begin();i!=tuple.end();i++)
        for(IntegerVectorList::const_iterator j=i->begin();j!=i->end();j++)
          {assert(d==j->size());}

      return d;
    }
    IntegerVectorList toNonSpecialSubspace(IntegerVectorList const &l, IntegerVector const &special)
    {
      int n=special.size();
      int nnonspecial=n-special.sum();
      IntegerVectorList ret;
      for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
        {
          IntegerVector temp(nnonspecial);
          int J=0;
          for(int j=0;j<n;j++)if(!special[j])temp[J++]=(*i)[j];
          ret.push_back(temp);
        }
      return ret;
    }
    void mainForProjection(IntegerVectorListList const &tuple, IntegerVector const *special=0)
    {
      int D=findD(tuple);
      int codim=coDimensionOfResultantVariety(tuple,D,0);

  //    if((special && (codim - special->sum()!=1))|| (!special && codim!=1))
      if(codim!=1)
        {
          debug<<"The --projection option only works for hypersurfaces.\n";
          assert(0);
        }

      IntegerMatrix A=cayleyConfiguration(tuple,findD(tuple));
      int N=A.getHeight();//big ambientdim
      int n=A.getHeight();
      if(special)n-=special->sum();
      vector<pair<int,int> > intervals=tupleToIntervals(tuple);

      IntegerVectorList nonSpecialGenerators;
      IntegerVectorList specialGenerators;
      if(special)
        for(int i=0;i<special->size();i++)
          if(!(*special)[i])
            nonSpecialGenerators.push_back(IntegerVector::standardVector(N,i));
          else
            specialGenerators.push_back(IntegerVector::standardVector(N,i));

      SelectionIterator iter(intervals);

      IntegerVectorList linealityGen=A.transposed().getRows();


      vector<PolyhedralCone> F;

      do
        {
          IntegerVectorList gen;
          int I=0;
          for(int i=0;i<iter.size();i++)for(int j=0;j<iter.sizeOfIth(i);j++,I++)if(!iter.chosen(i,j))gen.push_back(IntegerVector::standardVector(N,I));

          if(special)
            {
              IntegerVectorList temp=nonSpecialGenerators;
              for(IntegerVectorList::const_iterator i=linealityGen.begin();i!=linealityGen.end();i++)temp.push_back(*i);
              for(IntegerVectorList::const_iterator i=gen.begin();i!=gen.end();i++)temp.push_back(*i);

              if(::rank(rowsToIntegerMatrix(temp,N))!=N){log2 debug<<"Skipping\n";continue;}
            }


          PolyhedralCone C=PolyhedralCone::givenByRays(gen,linealityGen,N);
          if(special)
            {
              C=intersection(C,PolyhedralCone(IntegerVectorList(),specialGenerators,N));
              C=PolyhedralCone(toNonSpecialSubspace(C.getHalfSpaces(),*special),toNonSpecialSubspace(C.getEquations(),*special),n);
            }

/*          PolyhedralCone C=(special)?
              PolyhedralCone::givenByRays(toNonSpecialSubspace(gen,*special),toNonSpecialSubspace(linealityGen,*special),n)
            :
              PolyhedralCone::givenByRays(gen,linealityGen,n);*/
          C.canonicalize();
   //       debug<<toNonSpecialSubspace(gen,*special);
  //        debug<<toNonSpecialSubspace(linealityGen,*special);
          log2 debug<<":"<<C.dimension()<<"\n";
          if(C.dimension()==n-1)
            {
              log2 debug<<"adding\n";
              F.push_back(C);
            }
        }
      while(++iter);

      vector<BSPTree::SmallCone> a=BSPTree::buildPointers(F);
      BSPTree tree(n,a,0,false);
//      debug<<"Number of regions"<<tree.numberOfRegions()<<"\n";

      {
        BSPTreeTraverser2 traverser(n,tree,false);

        SymmetryGroup s(n);
        if(symmetryOption.getValue())
          {
            IntegerVectorList generators=FileParser(Stdin).parseIntegerVectorList();
            s.computeClosure(generators);
            s.createTrie();
          }
        SymmetricTargetFanBuilder target(n,s);
        symmetricTraverse(traverser,target,&s);
        target.getFanRef().printWithIndices(&pout,
            (symmetryOption.getValue()?FPF_group|FPF_conesCompressed:0)|
              (optionIgnoreCones.getValue()?0:FPF_conesExpanded)|
              FPF_maximalCones|FPF_cones/*|FPF_multiplicities*/,
              &s);
      }
    }
    Field f(){return Q;}
  int main()
  {
	  {//this code fails if we use shared pointers for ref counting:
	      	FieldElement a(Q.zHomomorphism(0));
	      	{
			  FieldVector A(a.getField(),0);
		  }
		  {
			  FieldVector B(a.getField(),0);
		  }
	  }
    IntegerVectorListList tuple;
    if(optionVectorInput.getValue())
      {
        tuple=FileParser(Stdin).parseIntegerVectorListList();
      }
    else
      {
        PolynomialSet g=FileParser(Stdin).parsePolynomialSetWithRing();

        tuple=g.exponents();
      }
    if(optionSpecialization.getValue() ||optionPrettyPrint.getValue())
      {
        if(optionPrettyPrint.getValue())
          {
            LatexPrinter P(Stdout);
            P.printVectorListList(tuple);
            P.printNewLine();
            return 0;
          }

        int D=findD(tuple);
        IntegerMatrix g;
        IntegerMatrix A=cayleyConfiguration(tuple,findD(tuple));//cayleyConfiguration(g);
        vector<pair<int,int> > intervals=tupleToIntervals(tuple);//polynomialSetToIntervals(g);
//        debug<<"Cayley configuration:"<<A.getRows()<<"\n";
        A=A.transposed();
        Triangulation2::makeConfigurationFullDimensional(A);
        A=A.transposed();
//         debug<<"Cayley configuration:"<<A.getRows()<<"\n";

        IntegerVector special=FileParser(Stdin).parseIntegerVector().supportAsZeroOneVector();
        if(optionProjection.getValue())
          {
            mainForProjection(tuple,&special);
          }
        else
          {
        IntegerVectorList exponentsL;
//        list<IntegerVectorList> exponentsA=g.exponents();

        IntegerVectorListList exponentsA=tuple;

        log1 debug<<"FRONT"<<exponentsA.front();
/*        vector<pair<int,int> > intervals;
        int I=0;
        for(list<IntegerVectorList>::const_iterator i=exponentsA.begin();i!=exponentsA.end())
          {
            intervals.push_back(pair<int,int>(I,I+i->size());
            I+=i->size();
          }*/
        for(list<IntegerVectorList>::const_iterator i=exponentsA.begin();i!=exponentsA.end();i++)for(IntegerVectorList::const_iterator j=i->begin();j!=i->end();j++)exponentsL.push_back(*j);

        IntegerMatrix exponentMatrix=rowsToIntegerMatrix(exponentsL).transposed();

        if(isSpecializedResultantEmpty(exponentMatrix,intervals,special))
          {
            debug<<"The specialized resultant variety is empty\n";
            assert(0);
          }

        /*
         * Find generic starting point omega.
         */

        PolyhedralCone toBeAvoided=intersection(specializedToSubspace(special),PolyhedralCone::givenByRays(IntegerVectorList(),A.transposed().getRows(),A.getHeight()));

        IntegerVectorList l=perturbationSequenceRek(exponentMatrix,intervals,IntegerVector::allOnes(exponentsL.size()),special,toBeAvoided);
        log1 D(l);
        IntegerMatrix cayley=cayleyConfiguration(exponentsA,D).transposed();
        Triangulation2::makeConfigurationFullDimensional(cayley);

        PolyhedralCone secondaryCone=perturbationSequenceToVectorInSecondaryCone(l, cayley);
        IntegerVector omega=intersection(secondaryCone,specializedToSubspace(special)).getRelativeInteriorPoint();

//        D(secondaryCone);
/*        {
          IntegerVector omega;
          int codimension=coDimensionOfResultantVariety(tuple,findD(tuple),&omega);
          if(secondaryCone.dimension()+codimension!=secondaryCone.ambientDimension())
        }
*/
        ResultantFanSpecializationTraverser traverser(tuple,D,A,intervals,special,omega,special.size()-special.sum()-coDimensionOfResultantVariety(tuple,findD(tuple),0));
        int n=traverser.refToPolyhedralCone().ambientDimension();
        SymmetryGroup s(n);

        SymmetricTargetFanBuilder target(n,s);
        symmetricTraverse(traverser,target,&s);
        target.getFanRef().printWithIndices(&pout,
            (symmetryOption.getValue()?FPF_group|FPF_conesCompressed:0)|
              (optionIgnoreCones.getValue()?0:FPF_conesExpanded)|
              FPF_maximalCones|FPF_cones,
              &s);
          }
        return 0;
      }



//    PolynomialSet g=FileParser(Stdin).parsePolynomialSetWithRing();

//    IntegerVectorListList tuple=g.exponents();
    log1 debug<<tuple;


//    assert(optionVectorInput.getValue()==false);
    if(!oldOption.getValue())
      {
        if(optionCodimension.getValue())
          {
            IntegerVector omega;
            int codimension=coDimensionOfResultantVariety(tuple,findD(tuple),&omega);

            pout<<codimension<<"\n";
            return 0;
          }
        else
        if(optionProjection.getValue())
          {
            mainForProjection(tuple);
          }
        else
          {
            IntegerMatrix A=cayleyConfiguration(tuple,findD(tuple));//cayleyConfiguration(g);
//            IntegerMatrix A=cayleyConfiguration(g);
            A=A.transposed();
            Triangulation2::makeConfigurationFullDimensional(A);
            A=A.transposed();
      //      debug<<"THECONF"<<A.getRows();
            //debug<<"Cayley configuration:"<<A.getRows()<<"\n";


      	  ResultantFanTraverser traverser(/*g*/tuple,A);
            int n=traverser.refToPolyhedralCone().ambientDimension();
            SymmetryGroup s(n);

            if(symmetryOption.getValue())
              {
                IntegerVectorList generators=FileParser(Stdin).parseIntegerVectorList();
                s.computeClosure(generators);
                s.createTrie();
              }

            SymmetricTargetFanBuilder target(n,s);
            symmetricTraverse(traverser,target,&s);
            target.getFanRef().printWithIndices(&pout,
                (symmetryOption.getValue()?FPF_group|FPF_conesCompressed:0)|
                  (optionIgnoreCones.getValue()?0:FPF_conesExpanded)|
                  FPF_maximalCones|FPF_cones|FPF_multiplicities,
                  &s);
          }
        return 0;
      }



    IntegerVector omega;
//    int codimension=coDimensionOfResultantVariety(g,&omega);
    int codimension=coDimensionOfResultantVariety(tuple,findD(tuple),&omega);

    log1 debug<<"Codimension of resultant variety:"<<codimension<<"\n";
    log1 debug<<omega<<"\n";


    IntegerMatrix A=cayleyConfiguration(tuple,findD(tuple));//cayleyConfiguration(g);
//    IntegerMatrix A=cayleyConfiguration(g).transposed();
    log1 debug<<"Cayley configuration:"<<A.getRows()<<"\n";
    Triangulation2::makeConfigurationFullDimensional(A);

    {
      IntegerVectorList M;
      for(int i=0;i<omega.size();i++)if(omega[i])M.push_back(IntegerVector::standardVector(omega.size(),i));
      MatrixTermOrder TO(M);
      Triangulation2 T(A);
      T.triangulate();
      T.changeToTriangulationInducedBy(TO);
      PolyhedralCone C=T.secondaryCone();
      C.canonicalize();
      IntegerVectorList inequalities=C.getHalfSpaces();
      IntegerVectorList equations=C.getEquations();
      for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
        if(coordinatewiseProduct(omega,*i).isZero())equations.push_back(*i);
      PolyhedralCone C2(inequalities,equations,C.ambientDimension());
      debug<<C2.getRelativeInteriorPoint();
    }


    int n=A.getWidth();

    SymmetryGroup s(n);

    Triangulation2 t(A);


    t.triangulate();

    SymmetricTargetFanBuilder target(n,s);

    {
      SecondaryFanTraverser traverser(t);
      symmetricTraverse(traverser,target,&s);
    }

    //////////////////////////////////////////////////////
    A=A.transposed();
    //////////////////////////////////////////////////////

    PolyhedralFan f2=target.getFanRef();

    for(int i=0;i<codimension;i++)
      f2=f2.facetComplex();

    PolyhedralFan f3(f2.getAmbientDimension());
    for(PolyhedralFan::coneIterator ii=f2.conesBegin();ii!=f2.conesEnd();ii++)
      {
      log1{static int t;t++;debug<<t<<"\n";}
        IntegerVector w=ii->getRelativeInteriorPoint();

        log1 {static int t;t++;if(!(t&127))debug<<t<<"\n";}

        {
          IntegerMatrix Atransposed=A.transposed();
          Triangulation2 t(Atransposed);

          /* Convert a Triangulation to a Triangulation2 */
          {
            list<Triangulation::Cone> T=Triangulation::triangulate(A);
            for(list<Triangulation::Cone>::const_iterator i=T.begin();i!=T.end();i++)
              {
                IntegerVector v=i->size();
                int J=0;
                for(Triangulation::Cone::const_iterator j=i->begin();j!=i->end();j++,J++)
                  v[J]=*j;
                  t.bases.insert(v);
              }
          }
          WeightReverseLexicographicTermOrder T(w);
          t.changeToTriangulationInducedBy(T);

//              t.print(debug);

          map<IntegerVector,set<int> > equivalenceClasses;
          for(set<IntegerVector>::const_iterator i=t.bases.begin();i!=t.bases.end();i++)
            {
              IntegerVectorList M;
              IntegerVector temp(1);
              for(int j=0;j<i->size();j++)
                {temp[0]=w[(*i)[j]];M.push_back(concatenation(A[(*i)[j]],temp));}
              IntegerMatrix M2=rowsToIntegerMatrix(M);
              IntegerVector normal=vectorInKernel(M2);
              for(int j=0;j<i->size();j++)equivalenceClasses[normal].insert((*i)[j]);
            }
          bool someOK=false;
          for(map<IntegerVector,set<int> >::const_iterator i=equivalenceClasses.begin();i!=equivalenceClasses.end();i++)
            {
              bool OK=true;
              int first=0;
              assert(0);//FIX THE CODE BELOW WHEN THERE IS NO g
/*              for(PolynomialSet::const_iterator j=g.begin();j!=g.end();j++)
                {
                  int last=first+j->numberOfTerms();
                  int numberFound=0;
                  for(int k=first;k<last;k++)if(i->second.count(k))numberFound++;
                  if(numberFound<2)OK=false;
                  first=last;
                }*/
              if(OK)someOK=true;
            }
          if(someOK)f3.insert(*ii);
        }
      }
    f3.printWithIndices(&pout);

    return 0;
  }
};

static ResultantFanApplication theApplication;

