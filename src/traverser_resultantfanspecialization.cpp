/*
 * traverser_resultantfanspecialization.cpp
 *
 *  Created on: Jul 1, 2011
 *      Author: anders
 */

#include "traverser_resultantfanspecialization.h"
#include "traverser_secondaryfan.h"
#include "halfopencone.h"




class Shrinking
{
public:
  IntegerVector subset;
  int newCount;
public:
  Shrinking(IntegerVector const &subset_):
    subset(subset_),
    newCount(subset_.sum())
  {
  }
  Shrinking(set<int> const &subsetIndices, int N):
    subset(N),
    newCount(subsetIndices.size())
    {
      for(set<int>::const_iterator i=subsetIndices.begin();i!=subsetIndices.end();i++)subset[*i]=1;
    }
  IntegerVector shrink(IntegerVector const &a)const
  {
    assert(a.size()==subset.size());
    IntegerVector ret(newCount);
    int I=0;
    for(int i=0;i<a.size();i++)if(subset[i])ret[I++]=a[i];
    return ret;
  }
  IntegerVectorList shrink(IntegerVectorList const &l)const
  {
    IntegerVectorList ret;
    for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
      ret.push_back(shrink(*i));
    return ret;
  }
  PolyhedralCone shrink(PolyhedralCone const &c)const
  {
    return PolyhedralCone(shrink(c.getHalfSpaces()),shrink(c.getEquations()),newCount);
  }
  IntegerVector expand(IntegerVector const &v)const
  {
    IntegerVector ret(subset.size());
    int I=0;
    for(int i=0;i<ret.size();i++)
      if(subset[i])ret[i]=v[I++];
    return ret;
  }
  IntegerVectorList expand(IntegerVectorList const &l)const
  {
    IntegerVectorList ret;
    for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
      ret.push_back(expand(*i));
    return ret;
  }
  // Adds a linear space
  PolyhedralCone expand(PolyhedralCone const &c)const
  {
    return PolyhedralCone(expand(c.getHalfSpaces()),expand(c.getEquations()),subset.size());
  }
  int newPosition(int old)
  {
//    if(old==subset.size())return newCount;
    int ret=0;
    old--;
    while(old>=0){if(subset[old])ret++;old--;}
    return ret;
  }
  int oldPosition(int n)
  {
    assert(n<newCount && n>=0);
    int ret=0;
    while(n>=0){if(subset[ret++])n--;}
    return ret-1;
  }
  vector<pair<int,int> > shrinkIntervals(vector<pair<int, int> >const &intervals)
    {
    vector<pair<int,int> > ret;
    for(vector<pair<int,int> >::const_iterator i=intervals.begin();i!=intervals.end();i++)
      ret.push_back(pair<int,int>(newPosition(i->first),newPosition(i->second)));
    return ret;
    }
};

bool isSpecializedResultantEmpty(IntegerMatrix const &exponents, vector<pair<int,int> > const &intervals, IntegerVector const &isSpecial)
{
  /*
   *  This check is done by creating a new point configuration whose resultant is full-dimensional iff the specialized resultant is non-empty.
   */
  list<list<IntegerVector> > extended;
  int nnonspecial=isSpecial.size()-isSpecial.sum();
  int I=0;
  int J=0;
  for(vector<pair<int,int> >::const_iterator i=intervals.begin();i!=intervals.end();i++)
    {
      list<IntegerVector> temp;
      for(int j=i->first;j!=i->second;j++,I++)
        {
          IntegerVector t(nnonspecial);
          if(!isSpecial[I])t[J++]=1;
          temp.push_back(concatenation(exponents.column(I),t));
        }
      extended.push_back(temp);
    }
  IntegerVector choice;
//    for(list<IntegerVectorList>::const_iterator i=extended.begin();i!=extended.end();i++)debug<<*i;
  return (coDimensionOfResultantVariety(extended,exponents.getHeight()+nnonspecial,&choice)!=0);
}


/**
 * The number of columns in the cayley configuration is n, which is also the number of entries of special and subset. Those are considered to be vectors of bools.
 * The subset tells the routine which columns of cayley to consider. The rest is simply ignored and zeroes are returned for those entries.
 * This routine will compute a vector in the resultant fan of the cayley, but outside the lineality space - that is the rowspace of cayley.
 */
IntegerVector nonTrivialVectorInSpecializedResultant(IntegerMatrix const &exponents, vector<pair<int,int> > const &intervals, IntegerVector const &subset, IntegerVector const &isSpecial)
{
  if(!(subset-IntegerVector::allOnes(subset.size())).isZero())
    {
      Shrinking s(subset);
      IntegerMatrix exponents2=exponents.submatrixColumnSubsetBoolean(subset);
      vector<pair<int,int> > intervals2=s.shrinkIntervals(intervals);
      IntegerVector subset2=IntegerVector::allOnes(subset.sum());
      IntegerVector isSpecial2=s.shrink(isSpecial);
      IntegerVector ret2=nonTrivialVectorInSpecializedResultant(exponents2,intervals2,subset2,isSpecial2);
      return s.expand(ret2);
    }
  // After this we may assume that subset is everything.

  assert(isSpecial.size()==exponents.getWidth());
  PolyhedralCone subspace=specializedToSubspace(isSpecial);

  // A will be the cayley embedding of the cofingurations - that is the matrix whose secodary fan we consider
  IntegerMatrix A=exponents;
  for(vector<pair<int,int> >::const_iterator i=intervals.begin();i!=intervals.end();i++)
    A.appendRow(IntegerVector::interval(subset.size(),i->first,i->second));

  log1 D(A.getRows());
  log1 D(isSpecial);
  // We compute the linality space the stable intersection
  PolyhedralCone linealityIntersection=intersection(PolyhedralCone::givenByRays(IntegerVectorList(),A.getRows(),A.getWidth()),subspace);
  int linealitySpaceDim=linealityIntersection.dimension();


  // Check if stable intersection is only the lineality space, by computing its dimension.
  // We need to transform the exponent matrix in to a list of lists of exponents,
  // then compute the dimension of this, and subtract the codimension of the subspace.
  {
    list<list<IntegerVector> > nonExtended;
    int I=0;
    int J=0;
    for(vector<pair<int,int> >::const_iterator i=intervals.begin();i!=intervals.end();i++)
      {
        list<IntegerVector> temp;
        for(int j=i->first;j!=i->second;j++,I++)
          temp.push_back(exponents.column(I));
        nonExtended.push_back(temp);
      }
    int bigCodimension=coDimensionOfResultantVariety(nonExtended,exponents.getHeight(),0);
    if(isSpecial.size()-(bigCodimension+isSpecial.sum())<=linealitySpaceDim)
      return IntegerVector(isSpecial.size());
  }


  // We now know that there is a vector outside the lineality space to be found.
  // We run through all edge choices.

  SelectionIterator iter(intervals);

  do
    {
      // For this edge choice, build the corresponding cone in the the tropical resultant
      IntegerVectorList gen;
      int I=0;
      for(int i=0;i<iter.size();i++)for(int j=0;j<iter.sizeOfIth(i);j++,I++)if(!iter.chosen(i,j))gen.push_back(IntegerVector::standardVector(exponents.getWidth(),I));
      PolyhedralCone C=PolyhedralCone::givenByRays(gen,A.getRows(),exponents.getWidth());
      C.canonicalize();

      // Compute sum and intersection with the subspace
      PolyhedralCone inter=intersection(C,subspace);
      PolyhedralCone su=sum(C,subspace);

      // Check if 1) there is a contribution to the stable intersection 2) the contribution
      // is more than the lineality space of the stable intersection
      if(su.dimension()==A.getWidth() && inter.dimension()>linealitySpaceDim)
        {
          // We find generators of the intersection cone (and linality space). One of these must be outside
          // the linality space of the stable intersection and can be returned.
          inter.canonicalize();
          IntegerVectorList candidates=inter.extremeRays();
          IntegerVectorList candidates2=inter.generatorsOfLinealitySpace();
          for(IntegerVectorList::const_iterator i=candidates2.begin();i!=candidates2.end();i++)candidates.push_back(*i);
          for(IntegerVectorList::const_iterator i=candidates.begin();i!=candidates.end();i++)
            {
              if(!linealityIntersection.contains(*i))
                return *i;
            }
          debug<<"NOT FOUND - THIS IS A BUG";
          assert(0);
        }
    }
  while(++iter);

  debug<<"NOT FOUND2 - THIS IS A BUG";
  assert(0);

  return IntegerVector(isSpecial.size());
}

/**
 * ToBeAvoided should be contained in the rowspace of the cayley + span of the standard vectors corresponding to
 * the variables not in the subproblem. It should also be contained in the non-special space.
 *
 * The to avoided space is a
 */
IntegerVectorList perturbationSequenceRek(IntegerMatrix const &exponents, vector<pair<int,int> > const &intervals, IntegerVector const &subset, IntegerVector const &isSpecial, PolyhedralCone const &toBeAvoided)
{
  log1 debug<<"INSUBSET"<<subset<<"\n";
  Shrinking s(subset);
  IntegerVector omega;

  //Compute cayley
/*  IntegerMatrix A=exponents.submatrixColumnSubsetBoolean(subset);
  {
    vector<pair<int,int> > intervals2=s.shrinkIntervals(intervals);
    for(vector<pair<int,int> >::const_iterator i=intervals2.begin();i!=intervals2.end();i++)
      A.appendRow(IntegerVector::interval(A.getWidth(),i->first,i->second));
  }
  */

  IntegerMatrix B=exponents;
  for(vector<pair<int,int> >::const_iterator i=intervals.begin();i!=intervals.end();i++)
    B.appendRow(IntegerVector::interval(B.getWidth(),i->first,i->second));
  for(int k=0;k<subset.size();k++)if(!subset[k])B.appendRow(IntegerVector::standardVector(subset.size(),k));
  PolyhedralCone B1=PolyhedralCone::givenByRays(IntegerVectorList(),B.getRows(),B.getWidth());
  PolyhedralCone intersectionOfSomething=intersection(B1,specializedToSubspace(isSpecial));

  if(toBeAvoided.dimension()==intersectionOfSomething.dimension())
    {
      omega=nonTrivialVectorInSpecializedResultant(exponents,intervals,subset,isSpecial);
    }
  else
    {
      log1 D(toBeAvoided);
      log1 D(intersectionOfSomething);
      IntegerVectorList s=intersectionOfSomething.generatorsOfLinealitySpace();
      for(IntegerVectorList::const_iterator i=s.begin();i!=s.end();i++)
        if(!toBeAvoided.contains(*i))
          {
            omega=*i;
            break;
          }
    }
  if(omega.isZero())
    {
//      debug<<"RETURNING EMPTY\n";
//      return IntegerVectorList();
      return toBeAvoided.generatorsOfLinealitySpace();
    }

  // We compute the subconfigurations which describe the link, but first we shrink the matrix of exponents (and
  // add the Cayley rows) since the Triangulation2 class does not allow columns to be ignored.
  IntegerMatrix newCayley=exponents.submatrixColumnSubsetBoolean(subset);
  vector<pair<int,int> > intervals2=s.shrinkIntervals(intervals);
  for(vector<pair<int,int> >::const_iterator i=intervals.begin();i!=intervals.end();i++)newCayley.appendRow(s.shrink(IntegerVector::interval(exponents.getWidth(),i->first,i->second)));

  Triangulation2 theTriangulation(newCayley);
  log1 debug<<"omega"<<omega<<"\n";

  IntegerVector newOmega=s.shrink(omega);
  WeightReverseLexicographicTermOrder T(newOmega);
  theTriangulation.triangulate();
  theTriangulation.changeToTriangulationInducedBy(T);
  PolyhedralCone secondaryCone=theTriangulation.secondaryCone().faceContaining(newOmega);

  PolyhedralCone newToBeAvoided=intersection(specializedToSubspace(isSpecial),s.expand(secondaryCone.span()));


  set<set<int> > subproblems=mixedCells(intervals2,theTriangulation,newCayley.transposed(),newOmega);

  log1 D((int)subproblems.size());


  // At this point we forgot the codimension of the resultant fan and need to recompute it.
  // Does this really handle specialization correctly?
  int minimalcodim=100;
  for(set<set<int> >::const_iterator i=subproblems.begin();i!=subproblems.end();i++)
    {
      IntegerVector newSubset(subset.size());
      for(set<int>::const_iterator j=i->begin();j!=i->end();j++)newSubset[s.oldPosition(*j)]=1;

      {
//          D(exponents.getRows());
//          D(newSubset);
          Shrinking S(newSubset);
          IntegerMatrix exponents2=rowsToIntegerMatrix(S.shrink(exponents.getRows()),S.newCount);
//          D(exponents2.getRows());

          IntegerVectorListList tuple2;
          vector<pair<int,int> >intervals2=S.shrinkIntervals(intervals);
          for(vector<pair<int,int> >::const_iterator i=intervals2.begin();i!=intervals2.end();i++)
            {
              IntegerVectorList temp;
              for(int j=i->first;j<i->second;j++)temp.push_back(exponents2.transposed()[j]);
              tuple2.push_back(temp);
            }
  //        D(tuple2);

          int codimension=coDimensionOfResultantVariety(tuple2,exponents.getHeight(),0);
//          D(codimension);
          if(codimension<minimalcodim)minimalcodim=codimension;
      }
    }
//  D(minimalcodim);

  // We now iterate through all sub configurations and pick one which contributes in the right dimension
  // Does this work correctly for specialised resultants?
  for(set<set<int> >::const_iterator i=subproblems.begin();i!=subproblems.end();i++)
    {
      IntegerVector newSubset(subset.size());
      for(set<int>::const_iterator j=i->begin();j!=i->end();j++)newSubset[s.oldPosition(*j)]=1;
      {
//          D(exponents.getRows());
//          D(newSubset);
          Shrinking S(newSubset);
          IntegerMatrix exponents2=rowsToIntegerMatrix(S.shrink(exponents.getRows()),S.newCount);

          IntegerVectorListList tuple2;
          vector<pair<int,int> >intervals2=S.shrinkIntervals(intervals);
          for(vector<pair<int,int> >::const_iterator i=intervals2.begin();i!=intervals2.end();i++)
            {
              IntegerVectorList temp;
              for(int j=i->first;j<i->second;j++)temp.push_back(exponents2.transposed()[j]);
              tuple2.push_back(temp);
            }
//          D(tuple2);

          int codimension=coDimensionOfResultantVariety(tuple2,exponents.getHeight(),0);
//          D(codimension);
          if(codimension!=minimalcodim)continue;
      }

      log1 debug<<"Iterating\n";
      IntegerVectorList ret=perturbationSequenceRek(exponents,intervals,newSubset,isSpecial,newToBeAvoided);
      // Perturb the returned vector to make it generic in the big resultant:
//      for(int k=0;k<exponents.getWidth();k++)if(subset[k])if(!isSpecial[k])if(i->count(k)==0)ret.push_back(IntegerVector::standardVector(subset.size(),k));
      // We add to our found omega epsilon times the vector found recursively and return:
      ret.push_front(omega);
//      debug<<"RETURNING"<<ret;
      return ret;
    }
  assert(0);
  return IntegerVectorList();
}


PolyhedralCone perturbationSequenceToVectorInSecondaryCone(IntegerVectorList const &l, IntegerMatrix const &configuration)
{
  MatrixTermOrder TO(l);
  Triangulation2 t(configuration);
  t.triangulate();
  t.changeToTriangulationInducedBy(TO);
  return t.secondaryCone().faceContainingPerturbed(l);
}


PolyhedralCone specializedToSubspace(IntegerVector const &isSpecial)
{
  IntegerVectorList temp;
  for(int i=0;i<isSpecial.size();i++)if(isSpecial[i])temp.push_back(IntegerVector::standardVector(isSpecial.size(),i));
  PolyhedralCone ret(IntegerVectorList(),temp,isSpecial.size());
  ret.canonicalize();
  return ret;
}


ResultantFanSpecializationTraverser::ResultantFanSpecializationTraverser(IntegerVectorListList const &tuple_,int D_, IntegerMatrix const &cayley, vector<pair<int,int> > intervals_, IntegerVector const &isSpecial_, IntegerVector omega, int dim):
		ConeTraverser(cayley.getHeight()),
  theTriangulation(cayley.transposed()),
  theCone(0),
  theConfiguration(cayley),
//  g(g_),
  tuple(tuple_),
  isSpecial(isSpecial_),
  intervals(intervals_),
  subspace(specializedToSubspace(isSpecial_)),
  D1(D_)
{
  log1 debug<<"OMEGA"<<omega;

  n=theTriangulation.getN();
    IntegerVector omega2;
//    int codimension=coDimensionOfResultantVariety(g,&omega2);
    int codimension=coDimensionOfResultantVariety(tuple,D1,&omega2);
    d=n-codimension-isSpecial.sum();
    IntegerVectorList M;
    M.push_back(omega);
    for(int i=0;i<omega.size();i++)if(omega[i])M.push_back(-IntegerVector::standardVector(omega.size(),i));
    MatrixTermOrder TO(M);
    theTriangulation.triangulate();
    theTriangulation.changeToTriangulationInducedBy(TO);
    PolyhedralCone C=theTriangulation.secondaryCone();
    C.canonicalize();
    assert(theCone.contains(omega));
    theCone=intersection(subspace,C.faceContaining(omega));
    theCone.canonicalize();

    assert(theCone.contains(omega));

    if(theCone.dimension()!=dim)
      {
        debug<<"Failed to produce starting cone of right dimension.\n";
        debug<<"Produced cone has dimension "<<theCone.dimension()<<"\n";
        debug<<"Expected "<<dim<<"\n";
        assert(0);
      }


    omega=theCone.getRelativeInteriorPoint();
  //  WeightReverseLexicographicTermOrder TO2(omega);//HERE
    WeightTermOrder TO2(omega);//HERE
    theTriangulation.changeToTriangulationInducedBy(TO2);
}


void ResultantFanSpecializationTraverser::changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector)
{
  IntegerVectorList M;
//  debug<<"CHANGING ridge"<<ridgeVector<<"ray"<<rayVector<<"\n";
  M.push_back(ridgeVector);//HERE
  M.push_back(rayVector);//HERE
for(int i=0;i<n;i++)M.push_back(-IntegerVector::standardVector(n,i));
  MatrixTermOrder T(M);
  theTriangulation.changeToTriangulationInducedBy(T);
  PolyhedralCone C=theTriangulation.secondaryCone();
  //debug>>C;
  C.canonicalize();
/*  IntegerVectorList inequalities=C.getHalfSpaces();
  IntegerVectorList equations=C.getEquations();
  for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
    if((dotLong(ridgeVector,*i)==0)&&(dotLong(rayVector,*i)==0))equations.push_back(*i);
  theCone=PolyhedralCone(inequalities,equations,C.ambientDimension());
 */
  IntegerVectorList temp;temp.push_back(ridgeVector);temp.push_back(rayVector);
  theCone=C.faceContainingPerturbed(temp);

  //debug>>theCone;
  theCone.canonicalize();
  assert(theCone.contains(ridgeVector));
  assert(theCone.link(ridgeVector).contains(rayVector));

  theCone=intersection(theCone,subspace);
  theCone.canonicalize();
  assert(theCone.dimension()==d);

  assert(theCone.contains(ridgeVector));
  assert(theCone.link(ridgeVector).contains(rayVector));
}


IntegerVectorList ResultantFanSpecializationTraverser::link(IntegerVector const &ridgeVector)
{
  // debug<<"Ridge"<<ridgeVector<<"\n";
  PolyhedralCone secondaryCone=theTriangulation.secondaryCone().faceContaining(ridgeVector);
  log1 debug << "Dim:"<<secondaryCone.dimension()<<" Expected dim"<<d-1+isSpecial.sum()<<"\n";

  set<set<int> > subproblems=mixedCells(intervals,theTriangulation,theConfiguration,ridgeVector);
  set<IntegerVector> ret;
//  debug<<"\nNPROBLEMS:"<<(int)subproblems.size()<<"\n";
  for(set<set<int> >::const_iterator i=subproblems.begin();i!=subproblems.end();i++)
    {
/*      debug<<"\nPROBLEM";
    for(set<int>::const_iterator j=i->begin();j!=i->end();j++)
      {
        debug<<*j;
      }
*/
    IntegerVector subset(n);
    for(set<int>::const_iterator j=i->begin();j!=i->end();j++)subset[*j]=1;
    Shrinking S(subset);
    vector<pair<int,int> > intervals2=S.shrinkIntervals(intervals);

    vector<int> subproblem;for(set<int>::const_iterator j=i->begin();j!=i->end();j++)subproblem.push_back(*j);

#if 1
    // Compute link in stable intersection by computing a restricted secondary fan and checking its rays
    {
        IntegerVectorList subMatrix;
        IntegerVector subIsSpecial(i->size());
        int J=0;
        for(set<int>::const_iterator j=i->begin();j!=i->end();j++,J++)
          {
            subMatrix.push_back(theConfiguration[*j]);
            subIsSpecial[J]=isSpecial[*j];
          }
        IntegerMatrix subA=rowsToIntegerMatrix(subMatrix);
        PolyhedralCone subSubSpace=specializedToSubspace(subIsSpecial);

//        debug<<"SUBA"<<subA.getRows();
//        if(isSpecializedResultantEmpty(subA.transposed().submatrix(3,0,5,subA.getHeight()),intervals2,subIsSpecial)){debug<<"SKIPPING\n";continue;}

        subA=subA.transposed();
        Triangulation2::makeConfigurationFullDimensional(subA);

        IntegerVectorList rays;

        int n=subA.getWidth();
        Triangulation2 t(subA);
        t.triangulate();
        //Check if secondary fan has too high dimension
        IntegerVectorList linealitySpaceGenerators;
        PolyhedralCone temp(IntegerVectorList(),subA.getRows(),subA.getWidth());
        IntegerVectorList small=temp.generatorsOfLinealitySpace();
        for(IntegerVectorList::const_iterator i=small.begin();i!=small.end();i++)linealitySpaceGenerators.push_back(S.expand(*i));
        PolyhedralCone secondaryLineality=PolyhedralCone(IntegerVectorList(),linealitySpaceGenerators,ridgeVector.size());
        PolyhedralCone linealityIntersection=intersection(secondaryLineality,subspace);
        int itsDimension=linealityIntersection.dimension();
        //assert(itsDimension==d || itsDimension==d-1);// it is possible that this should not be an assertion, but that the subproblem should just be ignored, since it cannot contribute anyway.
        if(itsDimension==d)
          {
            linealityIntersection.canonicalize();
            PolyhedralCone temp=intersection(secondaryCone.span(),subspace);temp.canonicalize();
            IntegerVectorList rays2=inducedLink(linealityIntersection,temp);
            for(IntegerVectorList::const_iterator k=rays2.begin();k!=rays2.end();k++)
              {
                //HERE WE NEED TO ADD A CHECK THAT THE LINEALITY SPACE ACTUALLY IS IN THE STABLE INTERSECTION
                //FIX THIS
                ret.insert(*k);
              }
          }
        else
          {
            SymmetryGroup s(n);
            SymmetricTargetFanBuilder target(n,s);
            SecondaryFanTraverser traverser(triangulationWithFullDimensionalIntersection(t,subSubSpace),subSubSpace);
            symmetricTraverse(traverser,target,&s);

            PolyhedralFan temp=target.getFanRef();
            rays=temp.getRaysInPrintingOrder(&s);

            {for(IntegerVectorList::const_iterator i=rays.begin();i!=rays.end();i++)
              {
                log1 debug<<"RAY"<<*i<<"\n";
                WeightReverseLexicographicTermOrder TO(*i);
                t.changeToTriangulationInducedBy(TO);
                set<set<int> > theMixedCells=mixedCells(intervals2,t,subA.transposed(),*i);

                log1 debug<<"TESTINGAAA"<<(int)theMixedCells.size()<<"\n";
                {for(set<set<int> >::const_iterator i=theMixedCells.begin();i!=theMixedCells.end();i++)
                  {
                	log1 debug<<"\nMIXEd";
                	for(set<int>::const_iterator j=i->begin();j!=i->end();j++)
                	{
                		log1 debug<<*j;
                	}
                	log1 debug<<"\n";
                  }
                }
                if(theMixedCells.size())
                  {
                    // we now know that the ray is in the resultantant of the link, but not necessarily in the
                    // specialized resultant. We need to check if this is the case by checking if the specialized
                    // resultants of the subconfigurations are non-empty.
                    bool isInSpecializedResultant=false;
                    for(set<set<int> >::const_iterator j=theMixedCells.begin();j!=theMixedCells.end();j++)
                      {
                        Shrinking S(*j,n);
                        IntegerMatrix subA3=subA.submatrixColumnSubsetBoolean(S.expand(IntegerVector::allOnes(j->size())));
                        subA3=subA3.submatrix(intervals.size(),0,subA3.getHeight(),subA3.getWidth());
                        log1 debug<<"subA3:"<<subA3.getRows()<<"\n";
                        vector<pair<int,int> > intervals3=S.shrinkIntervals(intervals2);
                        //debug<<"subA3:"<<subA3.getRows()<<"\n";
                        IntegerVector isSpecial3=S.shrink(subIsSpecial);
                        log1 debug<<"isSpecial3:"<<isSpecial3<<"\n";
                        for(vector<pair<int,int> > ::const_iterator k=intervals3.begin();k!=intervals3.end();k++)
                          {
                            log1 debug<<k->first<<k->second<<"\n";
                          }
                       if(!isSpecializedResultantEmpty(subA3,intervals3,isSpecial3)){isInSpecializedResultant=true;break;}
                      }
                    if(isInSpecializedResultant)ret.insert(S.expand(*i));
                  }
                log1 debug<<"TESTINGBBB\n";
              }
            }
    }
    }
#else
    // compute linke in stable intersection as a projection
//    vector<pair<int,int> > intervals=subsetToIntervals(g,*i);
    SelectionIterator iterator(intervals2);
    int numberOfIterations=0;
    do
        {
          IntegerVectorList generators;
          for(int j=0;j<n;j++)if(!(i->count(j)))generators.push_back(IntegerVector::standardVector(n,j));
          int I=0;
          for(int i=0;i<iterator.size();i++)for(int j=0;j<iterator.sizeOfIth(i);j++,I++)if(!(iterator.chosen(i,j)))generators.push_back(IntegerVector::standardVector(n,subproblem[I]));
//           debug<<generators;
//          debug<<theConfiguration.transposed().getRows();
          PolyhedralCone c=PolyhedralCone::givenByRays(generators,theConfiguration.transposed().getRows(),n);
//          debug<<"Dimension of resultant cone"<<c.dimension()<<"\n";
//          debug<<"Dimension of subspace intersection"<<intersection(subspace,c).dimension()<<"\n";
          subspace.canonicalize();
          c.canonicalize();
          PolyhedralCone s=sum(subspace,c);
//          debug<<"Dimension of sum"<<s.dimension()<<"\n";

          if((s.dimension()==n) && (intersection(subspace,c).dimension()==d))
            {
              PolyhedralCone A=intersection(subspace,c).link(ridgeVector);A.canonicalize();
              PolyhedralCone B=theCone.faceContaining(ridgeVector).span();B.canonicalize();
              IntegerVectorList temp=inducedLink(A,B);
              for(IntegerVectorList::const_iterator i=temp.begin();i!=temp.end();i++)ret.insert(*i);
            }
          numberOfIterations++;
        }
      while(++iterator);
      debug<<"Number of iterations\n"<<numberOfIterations;
#endif

    }
  IntegerVectorList ret2;
  for(set<IntegerVector>::const_iterator i=ret.begin();i!=ret.end();i++)
    ret2.push_back(*i);
  // debug<<"LINK AT"<<ridgeVector<<":\n";
  // debug<<ret2;
  //pout<<"LINK AT"<<ridgeVector<<":\n";
  //pout<<ret2;
  return ret2;
//  intersection(secondaryCone,subspace).printAsFan(&debug);
//  assert(0);
}


PolyhedralCone &ResultantFanSpecializationTraverser::refToPolyhedralCone()
{
  return theCone;
}
