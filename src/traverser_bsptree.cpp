/*
 * traverser_bsptree.cpp
 *
 *  Created on: Aug 23, 2011
 *      Author: anders
 */

#include "traverser_bsptree.h"
#include "parser.h"

BSPTreeTraverser::BSPTreeTraverser(int n_, BSPTree const & tree_):
  ConeTraverser(n_),
  tree(tree_)
{
  IntegerVectorList omega=IntegerMatrix::identity(n).getRows();
  theCone=tree.region(omega);
  theCone.canonicalize();
}

void BSPTreeTraverser::changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector)
{
  IntegerVectorList omega;
  omega.push_back(ridgeVector);
  omega.push_back(rayVector);
//debug<<"OMEGA"<<omega;

  theCone=tree.region(omega);
            theCone.canonicalize();
//            assert(theCone.dimensionOfLinealitySpace()==5);
/*            {

              PolyhedralFan temp(n);
              temp.insert(theCone);
              temp.printWithIndices(&debug);
            }*/
//            debug<<theCone;
}


IntegerVectorList BSPTreeTraverser::link(IntegerVector const &ridgeVector)
{
  IntegerVector v=theCone.link(ridgeVector).getUniquePoint();
  IntegerVectorList ret;
  ret.push_back(v);
  ret.push_back(-v);
  return ret;
}

PolyhedralCone & BSPTreeTraverser::refToPolyhedralCone()
{
  return theCone;
}

bool BSPTreeTraverser::hasNoState()const
{
  return true;
}



BSPTreeTraverser2::BSPTreeTraverser2(int n_, BSPTree const & tree_, bool vertexConstructing_):
		ConeTraverser(n_),
  tree(tree_)
{
		vertexConstructing=vertexConstructing_;
  IntegerVector omega=tree.randomVectorInComplement();
//  IntegerVector omega=StringParser("(-132200,-80300,309800,-97300,462600,-301899,-320001,159300,138900,-60401,-78499)").parseIntegerVector();
  theCone=tree.regionFast(omega);
//  debug<<theCone;
//  assert(0);
  theCone.canonicalize();
}

void BSPTreeTraverser2::changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector)
{
//  D(ridgeVector);
//  D(rayVector);
  IntegerVector omega=tree.perturbationToVectorInComplement(ridgeVector,rayVector);


  PolyhedralCone old=theCone;

  theCone=tree.regionFast(omega);


  theCone.canonicalize();
/*
  if((old.faceContaining(ridgeVector)!=theCone.faceContaining(ridgeVector)))
    {
      D(ridgeVector);
      D(rayVector);
      D(old);
      debug>>old;
      D(theCone);
      debug>>theCone;
      assert(0);
    }*/
  IntegerVector lastFacetVector=ridgeVector;
  IntegerVectorList temp=old.faceContaining(ridgeVector).getEquations();
  assert(temp.size()==1);
  IntegerVector lastFlipDirection=temp.front();

  if(dotLong(lastFlipDirection,rayVector)<0)lastFlipDirection=-lastFlipDirection;

//  debug<<tree.multiplicity(lastFacetVector,lastFlipDirection) <<" "<<lastFlipDirection<<"\n";

  currentVertexCoordinates+=tree.multiplicity(lastFacetVector,lastFlipDirection)*lastFlipDirection;
}


IntegerVectorList BSPTreeTraverser2::link(IntegerVector const &ridgeVector)
{
  IntegerVector v=theCone.link(ridgeVector).getUniquePoint();
  IntegerVectorList ret;
  ret.push_back(v);
  ret.push_back(-v);
  return ret;
}

PolyhedralCone & BSPTreeTraverser2::refToPolyhedralCone()
{
  return theCone;
}

bool BSPTreeTraverser2::hasNoState()const
{
	return false;
//  return true;
}

