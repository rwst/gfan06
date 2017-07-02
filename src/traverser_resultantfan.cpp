/*
 * traverser_resultantfan.cpp
 *
 *  Created on: Apr 29, 2011
 *      Author: anders
 */

#include "traverser_resultantfan.h"
#include "traverser_secondaryfan.h"
#include <iostream>
#include "halfopencone.h"
#include "regularsubdivision.h"
#include "log.h"

IntegerMatrix cayleyConfiguration(PolynomialSet const &g)
{
  int ambientDim=g.getRing().getNumberOfVariables();
  int numberOfPolytopes=g.size();
  IntegerVectorList L;
  int I=0;
  for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++,I++)
    {
    for(TermMap::const_iterator j=i->terms.begin();j!=i->terms.end();j++)
      {
        L.push_back(concatenation(IntegerVector::standardVector(numberOfPolytopes,I),j->first.exponent));
      }
    }
    return rowsToIntegerMatrix(L);
}

IntegerMatrix cayleyConfiguration(list<IntegerVectorList> const &g, int d)
{
  int ambientDim=d;
  int numberOfPolytopes=g.size();
  IntegerVectorList L;
  int I=0;
  for(list<IntegerVectorList>::const_iterator i=g.begin();i!=g.end();i++,I++)
    {
    for(IntegerVectorList::const_iterator j=i->begin();j!=i->end();j++)
      {
        L.push_back(concatenation(IntegerVector::standardVector(numberOfPolytopes,I),*j));
      }
    }
  return rowsToIntegerMatrix(L);
}



set<set<int> > mixedCells(vector<pair<int,int> > const&intervals/*PolynomialSet const &g*/, Triangulation2 const &t, IntegerMatrix const &A, IntegerVector const &w)
{
  set<set<int> > ret;
#if 1
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
  for(map<IntegerVector,set<int> >::const_iterator i=equivalenceClasses.begin();i!=equivalenceClasses.end();i++)
#else
  set<set<int> > temp=regularSubdivision(A,w);
  for(set<set<int> >::const_iterator i=temp.begin();i!=temp.end();i++)

#endif
    {
      set<int> allPoints;
      for(int j=0;j<A.getHeight();j++)
        {
          IntegerVector temp(1);
          temp[0]=w[j];
          if(dotLong(i->first,concatenation(A[j],temp))==0)allPoints.insert(j);
        }
//      for(set<int>::const_iterator j=i->second.begin();j!=i->second.end();j++)debug<<*j<<",";
//      debug<<"\n";
      bool OK=true;
/*      int first=0;
      for(PolynomialSet::const_iterator j=g.begin();j!=g.end();j++)
      {
        int last=first+j->numberOfTerms();
*/
      for(int j=0;j<intervals.size();j++)
        {
          int first=intervals[j].first;
          int last=intervals[j].second;
//          debug<<first<<"-"<<last<<"\n";
        int numberFound=0;
//        for(int k=first;k<last;k++)if(i->second.count(k))numberFound++;
        for(int k=first;k<last;k++)if(allPoints.count(k))numberFound++;
        if(numberFound<2)OK=false;
        first=last;
      }
      if(OK)
        {
//        ret.insert(*i->second);
        ret.insert(allPoints);
        }
    }
  return ret;
}


vector<pair<int,int> > polynomialSetToIntervals(PolynomialSet const &g)
{
  vector<pair<int,int> > intervals;
  int index=0;
  for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    {
      intervals.push_back(pair<int,int>(index,index+i->numberOfTerms()));
      index+=i->numberOfTerms();
    }
  return intervals;
}


vector<pair<int,int> > tupleToIntervals(IntegerVectorListList const &tuple)
{
  vector<pair<int,int> > intervals;
  int index=0;
  for(IntegerVectorListList::const_iterator i=tuple.begin();i!=tuple.end();i++)
    {
      intervals.push_back(pair<int,int>(index,index+i->size()));
      index+=i->size();
    }
  return intervals;
}


set<set<int> > mixedCells(PolynomialSet const &g, Triangulation2 const &t, IntegerMatrix const &A, IntegerVector const &w)
{
  return mixedCells(polynomialSetToIntervals(g),t,A,w);
}

set<set<int> > mixedCells(IntegerVectorListList const &tuple, Triangulation2 const &t, IntegerMatrix const &A, IntegerVector const &w)
{
  return mixedCells(tupleToIntervals(tuple),t,A,w);
}

static void printCone(PolyhedralCone c)
{
  c.canonicalize();
  PolyhedralFan f(c.ambientDimension());
  f.insert(c);
  f.printWithIndices(&debug);
}


// contributed by Josephine Yu
int weightVectorToMultiplicity(/*PolynomialSet const &g,*/IntegerVectorListList const &tuple, IntegerMatrix const &theConfiguration, Triangulation2 const &theTriangulation, IntegerVector const &omega)
{
  int mult = 0;
  set<set<int> > theMixedCells = mixedCells(tuple,theTriangulation,theConfiguration,omega);
  for(set<set<int> >::const_iterator i=theMixedCells.begin();i!=theMixedCells.end();i++)
    {
      // should assert that each set contains exactly two elements

                //pick out submatrix whose ROWS are the vertices of the Cayley configuration in the mixed cell
                IntegerVectorList submatrix;
                for(set<int>::const_iterator j=i->begin();j!=i->end();j++){submatrix.push_back(theConfiguration[*j]);}
                IntegerMatrix newConfiguration=rowsToIntegerMatrix(submatrix,theConfiguration.getWidth());
                FieldMatrix hermiteForm=integerMatrixToFieldMatrix(newConfiguration, Q);
                toInteger(hermiteForm[0][0])+toInteger(hermiteForm[1][1]);
                hermiteForm.reduce(false,true);
                int latticeIndex =1;
                for(int k=0; k < hermiteForm.getWidth(); k++) {
                latticeIndex *= toInteger(hermiteForm[k][k]);}
                mult += abs(latticeIndex);
        }
        return mult;
}

static int findD(IntegerVectorListList const &tuple)
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

ResultantFanTraverser::ResultantFanTraverser(/*PolynomialSet const &g_,*/ IntegerVectorListList const &tuple_, IntegerMatrix const &theConfiguration_):
		ConeTraverser(theConfiguration_.getHeight()),
		  theTriangulation(theConfiguration_.transposed()),
  theCone(0),
//  g(g_),
  tuple(tuple_),
  theConfiguration(theConfiguration_)
{
	assert(theTriangulation.getN()==theConfiguration_.getHeight());
  assert(::rank(theConfiguration)==theConfiguration.getWidth());
  n=theTriangulation.getN();
  IntegerVector omega;
  int codimension=coDimensionOfResultantVariety(tuple,findD(tuple),&omega);
//  int codimension=coDimensionOfResultantVariety(g,&omega);
  d=n-codimension;

  IntegerVectorList M;
  M.push_back(omega);
  for(int i=0;i<omega.size();i++)if(omega[i])M.push_back(-IntegerVector::standardVector(omega.size(),i));
  MatrixTermOrder TO(M);
  theTriangulation.triangulate();
  theTriangulation.changeToTriangulationInducedBy(TO);
//  debug<<"omega"<<omega<<"\n";
//  debug<<"mult:"<<weightVectorToMultiplicity(tuple,theConfiguration,theTriangulation,omega)<<"\n";
  PolyhedralCone C=theTriangulation.secondaryCone();
  C.canonicalize();
//  debug<<"mult:"<<weightVectorToMultiplicity(tuple,theConfiguration,theTriangulation,omega)<<"\n";
  assert(theCone.contains(omega));
  IntegerVectorList inequalities=C.getHalfSpaces();
  IntegerVectorList equations=C.getEquations();
  for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
    if(coordinatewiseProduct(omega,*i).isZero())equations.push_back(*i);
  theCone=PolyhedralCone(inequalities,equations,C.ambientDimension());
  theCone.canonicalize();

  assert(theCone.contains(omega));
  theCone.setMultiplicity(weightVectorToMultiplicity(tuple,theConfiguration,theTriangulation,omega));

  omega=theCone.getRelativeInteriorPoint();
//  WeightReverseLexicographicTermOrder TO2(omega);//HERE
  WeightTermOrder TO2(omega);//HERE
  theTriangulation.changeToTriangulationInducedBy(TO2);

//  printCone(theCone);
}



void ResultantFanTraverser::changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector)
{
  IntegerVectorList M;
//  debug<<"CHANGING ridge"<<ridgeVector<<"ray"<<rayVector<<"\n";
  M.push_back(ridgeVector);//HERE
  M.push_back(rayVector);//HERE
for(int i=0;i<n;i++)M.push_back(-IntegerVector::standardVector(n,i));
  MatrixTermOrder T(M);
  theTriangulation.changeToTriangulationInducedBy(T);
  PolyhedralCone C=theTriangulation.secondaryCone();
  C.canonicalize();
  IntegerVectorList inequalities=C.getHalfSpaces();
  IntegerVectorList equations=C.getEquations();
  for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
    if((dotLong(ridgeVector,*i)==0)&&(dotLong(rayVector,*i)==0))equations.push_back(*i);
  theCone=PolyhedralCone(inequalities,equations,C.ambientDimension());
  theCone.canonicalize();
  theCone.setMultiplicity(weightVectorToMultiplicity(tuple,theConfiguration,theTriangulation,theCone.getRelativeInteriorPoint()));

  assert(theCone.contains(ridgeVector));
  assert(theCone.link(ridgeVector).contains(rayVector));
}

vector<pair<int,int> > subsetToIntervals(PolynomialSet const &g, set<int> const &subset)
    {//improve complexity
  vector<pair<int,int> > intervals;
  int index=0;
  int indexOriginal=0;
  for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    {
      int nt=0;
      for(int j=indexOriginal;j<indexOriginal+i->numberOfTerms();j++)nt+=subset.count(j);
      intervals.push_back(pair<int,int>(index,index+nt));
      index+=nt;
      indexOriginal+=i->numberOfTerms();
    }
  return intervals;
    }
vector<pair<int,int> > subsetToIntervals(IntegerVectorListList const &tuple, set<int> const &subset)
    {//improve complexity
  vector<pair<int,int> > intervals;
  int index=0;
  int indexOriginal=0;
  for(IntegerVectorListList::const_iterator i=tuple.begin();i!=tuple.end();i++)
    {
      int nt=0;
      for(int j=indexOriginal;j<indexOriginal+i->size();j++)nt+=subset.count(j);
      intervals.push_back(pair<int,int>(index,index+nt));
      index+=nt;
      indexOriginal+=i->size();
    }
  return intervals;
    }

IntegerVectorList quickLink(IntegerMatrix const &m, set<int> const &cell, vector<pair<int,int> > const &intervals, IntegerVectorList const &linealitySpace)
{
  int n=m.getHeight();
  IntegerVectorList ret;
  set<int> zeroHeight;
  set<int> otherHeight;

  for(vector<pair<int,int> >::const_iterator i=intervals.begin();i!=intervals.end();i++)
    {
      set<int> inInterval;
      for(set<int>::const_iterator j=cell.begin();j!=cell.end();j++)
        if(*j>=i->first && *j<i->second)inInterval.insert(*j);
      if(inInterval.size()==2)
        for(set<int>::const_iterator j=inInterval.begin();j!=inInterval.end();j++)zeroHeight.insert(*j);
      else
        for(set<int>::const_iterator j=inInterval.begin();j!=inInterval.end();j++)otherHeight.insert(*j);
    }
  assert(otherHeight.size()==3);

  for(set<int>::const_iterator i=otherHeight.begin();i!=otherHeight.end();i++)
    {
/*      IntegerVector v(n)=IntegerVector::all(n);

      for(set<int>::const_iterator j=otherHeight.begin();j!=otherHeight.end();j++)
        if(j!=i)
          v[*j]=1;
*/      IntegerVectorList rays;
//  for(int j=0;j<n;j++)if(!zeroHeight.count(j))if(!otherHeight.count(j))rays.push_back(IntegerVector::standardVector(n,j));
      rays.push_back(IntegerVector::standardVector(n,*i));
//      IntegerVectorList lines=m.transposed().getRows();
      PolyhedralCone c=PolyhedralCone::givenByRays(rays,linealitySpace,n);
      c.canonicalize();
      IntegerVector w=c.getUniquePoint();
      if(!w.isZero())ret.push_back(w);
    }
  return ret;
}

IntegerVectorList ResultantFanTraverser::link(IntegerVector const &ridgeVector)
{
//  debug<<"THECONE.int"<<theCone.getRelativeInteriorPoint();
//  debug<<"RIDGE"<<ridgeVector<<"\n";

  IntegerVectorList linealitySpace=theCone.faceContaining(ridgeVector).generatorsOfLinealitySpace();

  set<IntegerVector> ret;
  set<set<int> > subproblems=mixedCells(tuple,theTriangulation,theConfiguration,ridgeVector);//HERE
  assert(!subproblems.empty());
  for(set<set<int> >::const_iterator i=subproblems.begin();i!=subproblems.end();i++)
    {
    //pick out submatrix
    IntegerVectorList submatrix;
    for(set<int>::const_iterator j=i->begin();j!=i->end();j++)submatrix.push_back(theConfiguration[*j]);
    IntegerMatrix newConfiguration=rowsToIntegerMatrix(submatrix,theConfiguration.getWidth()).transposed();

//    debug<<newConfiguration.getRows();

    //update intervals
    vector<pair<int,int> > intervals=subsetToIntervals(/*g*/tuple,*i);

    //compute its secondary fan
    Triangulation2::makeConfigurationFullDimensional(newConfiguration);

    if(newConfiguration.getHeight()==newConfiguration.getWidth()-(n-d))
      {
      //debug<<"SPECIAL\n";
      //We are in the case where the lineality space of the secondary fan of the subproblem has
      //codimension equal to that of the resultant variety. In this case the contribution of the
      //subproblem is two rays, which


      IntegerMatrix M(theConfiguration.getHeight(),theConfiguration.getWidth());
      for(set<int>::const_iterator j=i->begin();j!=i->end();j++)M[*j]=theConfiguration[*j];

      IntegerVectorList temp=M.transposed().getRows();
      for(int j=0;j<n;j++)if(!i->count(j))temp.push_back(IntegerVector::standardVector(n,j));
      PolyhedralCone span=PolyhedralCone::givenByRays(IntegerVectorList(),temp,n);
//      PolyhedralCone span=PolyhedralCone(IntegerVectorList(),M.transposed().getRows(),theCone.ambientDimension());
//      PolyhedralCone span=PolyhedralCone(IntegerVectorList(),theCone.getEquations(),theCone.ambientDimension());
      PolyhedralCone C=theCone.link(ridgeVector).linealitySpace().dualCone();
      PolyhedralCone D=intersection(span,C);
      assert(D.dimension()==1);
      IntegerVector v=normalized(D.generatorsOfLinealitySpace().front());
      ret.insert(v);
      ret.insert(-v);
      continue;
      }
       // debug<<"NONSPECIAL\n";
#if 0
//    debug<<"NONSPECIAL\n";

//    debug<<newConfiguration.getRows();
    Triangulation2 T2(newConfiguration);
    T2.triangulate();
    SecondaryFanTraverser traverser(T2);
    SymmetryGroup trivial(submatrix.size());
    SymmetricTargetFanBuilder target(submatrix.size(),trivial);
    symmetricTraverse(traverser,target);

    newConfiguration=newConfiguration.transposed();

//    debug<<newConfiguration.getRows();

    //for each ray check if it has a mixed cell
//    bool noneAdded=true;
    IntegerVectorList toBeLifted;
    IntegerVectorList rays=target.getFanRef().getRaysInPrintingOrder(0,false);
    for(IntegerVectorList::const_iterator j=rays.begin();j!=rays.end();j++)
      {
        WeightReverseLexicographicTermOrder T(*j);
        T2.changeToTriangulationInducedBy(T);
        set<set<int> > cells2=mixedCells(intervals,T2,newConfiguration,*j);

        //if so lift the vector back into full space and add to ret
        if(!cells2.empty())toBeLifted.push_back(*j);
      }
#else
//    debug<<"To be lifted:"<<toBeLifted;
    set<int> temp;for(int k=0;k<i->size();k++)temp.insert(k);
    IntegerVectorList toBeLifted=quickLink(newConfiguration.transposed(),temp,subsetToIntervals(/*g*/tuple,*i),newConfiguration.getRows());
#endif
    for(IntegerVectorList::const_iterator j=toBeLifted.begin();j!=toBeLifted.end();j++)
      {
        IntegerVector lift(n);
        int K=0;
        for(set<int>::const_iterator k=i->begin();k!=i->end();k++,K++)lift[*k]=(*j)[K];
        ret.insert(lift);
//        debug<<lift<<"\n";
//        noneAdded=false;
      }

    /*    if(noneAdded)
      {
        //WE NEED TO FIGURE OUT WHAT THE REASON FOR THIS HAPPENING IS
        IntegerVector v=theCone.link(ridgeVector).getUniquePoint();
        ret.insert(v);
        ret.insert(-v);
      }*/
    }

  IntegerVectorList ret2;
  for(set<IntegerVector>::const_iterator i=ret.begin();i!=ret.end();i++)ret2.push_back(*i);

//  debug<<"LINK AT "<<ridgeVector<<":\n"<<ret2;
  assert(!ret2.empty());
  return ret2;
}

PolyhedralCone & ResultantFanTraverser::refToPolyhedralCone()
{
  return theCone;
}
