/*
 * bsptree.h
 *
 *  Created on: Feb 4, 2011
 *      Author: anders
 */

#ifndef BSPTREE_H_
#define BSPTREE_H_


#include <vector>
#include <map>

#include "polyhedralcone.h"
#include "symmetry.h"

/**
 * The BSPTree is binary space partition tree. Its main feature is that it can be used to recover convex connected
 * components of complements to unions of cones. However, with the algorithm of [Jensen, Yu], it is no longer an
 * advantage to create the tree. In fact, doing so in higher dimensions is bad. Therefore the tree is not constructed
 * until buildTree is called. Notice that this can happen in the constructor.
 *
 * The main functionality of the BSPTree class is therefore now to store symmetric fans compactly via the SmallCone class.
 * The method regionFast() may be called even when the tree has not been build.
 */
class BSPTree
{
public:
#if 1
    typedef struct SmallCone{
    IntegerVector permutation;
    PolyhedralCone const *theCone;
    SmallCone(PolyhedralCone const *theCone_, IntegerVector const &permutation_):
      permutation(permutation_),
      theCone(theCone_)
    {
    }
    IntegerVectorList getEquations()const
    {
      IntegerVectorList ret;SymmetryGroup::appendPermutedIntegerVectorList(theCone->getEquations(),permutation,ret);return ret;
      //return SymmetryGroup::permuteIntegerVectorList(theCone->getEquations(),permutation);
    }
    void assignEquation(IntegerVector &dest)const
    {
      IntegerVectorList const &l=theCone->getEquations();
      assert(l.size()==1);
      permutation.composeAssign(l.front(),dest);
//      SymmetryGroup::composeAssign(permutation,l.front(),dest);
    }
    IntegerVectorList getHalfSpaces()const
    {
      IntegerVectorList ret;SymmetryGroup::appendPermutedIntegerVectorList(theCone->getHalfSpaces(),permutation,ret);return ret;
//      return SymmetryGroup::permuteIntegerVectorList(theCone->getHalfSpaces(),permutation);
    }
    IntegerVectorList generatorsOfLinealitySpace()const
    {
      IntegerVectorList ret;SymmetryGroup::appendPermutedIntegerVectorList(theCone->generatorsOfLinealitySpace(),permutation,ret);return ret;
//      return SymmetryGroup::permuteIntegerVectorList(theCone->generatorsOfLinealitySpace(),permutation);
    }
    void appendGeneratorsOfLinealitySpace(IntegerVectorList &ret)const
    {
      SymmetryGroup::appendPermutedIntegerVectorList(theCone->generatorsOfLinealitySpace(),permutation,ret);
    }
    IntegerVectorList extremeRays()const
    {
      IntegerVectorList ret;SymmetryGroup::appendPermutedIntegerVectorList(theCone->extremeRays(),permutation,ret);return ret;
//      return SymmetryGroup::permuteIntegerVectorList(theCone->extremeRays(),permutation);
    }
    void appendExtremeRays(IntegerVectorList &ret)const
    {
      SymmetryGroup::appendPermutedIntegerVectorList(theCone->extremeRays(),permutation,ret);
    }
    bool contains(IntegerVector const &v)const
    {
      return theCone->contains(SymmetryGroup::composeInverse(permutation,v));
    }
    bool containsPerturbed(IntegerVectorList const &l)const
    {
    	if(!contains(l.front()))return false;
    	return theCone->containsPerturbed(SymmetryGroup::permuteInverseIntegerVectorList(permutation,l));
    }
    bool doesSatisfyInequalityExpensive(IntegerVector const &ineq)const
    {
      //IntegerVector ineq2=SymmetryGroup::composeInverse(permutation,ineq);
      static IntegerVector ineq2;
      permutation.composeInverseAssign(ineq,ineq2);
      return theCone->doesSatisfyInequalityExpensive(ineq2);
    }


    /*
     * 1: before
     * 0: same
     * -1: after
     */
    static bool beforeNonPert(IntegerVector const &u, IntegerVector const &p, IntegerVector const &h1, IntegerVector const &h2)
    {
      /*
       * The ray is (after scaling) parameterized as follows
       * v(t)=tu+b
       * where t goes from infinity through 0 to minus infinity.
       *
       * The times for intersection are
       * t_i=(hi*b)/(-u*hi)
       * WLOG we assume that u*hi>0
       *
       * Hence we check if
       * (h1*b)*(u*h2) < (h2*b)*(u*h1)
       */
      int64 uh1,uh2,ph1,ph2;
      dotLong4(u,p,h1,h2,uh1,uh2,ph1,ph2);
      bool flag=(uh1<0)^(uh2<0);
      if(ph1*uh2>ph2*uh1)return false^flag;
      if(ph1*uh2<ph2*uh1)return true^flag;
      return false;
    }

    /**
     * Check if line starting at u, passing through p intersects the cone in its relative interior.
     * Assumptions:
     * u is not in the hyperplane spanned by the cone.
     *
     * Return values
     * 2: passes through relative interior
     * 1: passes through boundary
     * 0: does not intersect
     */
    int doIntersectNonPerturbed(IntegerVector const &u, IntegerVector const &p)const
    {
//      IntegerVector u2=SymmetryGroup::composeInverse(permutation,u);
//      IntegerVector p2=SymmetryGroup::composeInverse(permutation,p);
      static IntegerVector u2;
      permutation.composeInverseAssign(u,u2);
      static IntegerVector p2;
      permutation.composeInverseAssign(p,p2);

      IntegerVector const &equation=theCone->getEquations().front();

      IntegerVectorList const &inequalities=theCone->getHalfSpaces(); //Does this copy the list?

      bool passesThroughClosed=true;
      bool passesThroughInterior=true;
      for(IntegerVectorList::const_iterator j=inequalities.begin();j!=inequalities.end();j++)
        {
          int64 ju=dotLong(*j,u2);

          if(ju==0)
            {
              int64 jp=dotLong(*j,p2);
              if(jp<=0)passesThroughInterior=false;
              if(jp<0)passesThroughClosed=false;
            }
          else if(ju>0)
            {
              int b=beforeNonPert(u2,p2,*j,equation);
              if(b==1)
                {
                  passesThroughInterior=false;
                  passesThroughClosed=false;
                }
              if(b==0)
                passesThroughInterior=false;
            }
          else
            {
              int b=beforeNonPert(u2,p2,equation,*j);
              if(b==1)
                 {
                   passesThroughInterior=false;
                   passesThroughClosed=false;
                 }
               if(b==0)
                 passesThroughInterior=false;
            }
          if(!passesThroughClosed && !passesThroughInterior)break;
        }
      if(passesThroughInterior)return 2;
      return (passesThroughClosed);
    }
    int multiplicity()const
    {
        return theCone->getMultiplicity();
    }
  }SmallCone;
#else
  typedef PolyhedralCone SmallCone;
#endif
  int n;
//  vector<PolyhedralCone> const &theCones;
  vector<SmallCone> const &theCones;
  vector<IntegerVector> hyperplanes;
  class Node
  {
  public:
    Node *leftRight[2];
    IntegerVector normal;
    vector<int> conesInHyperplane;
    Node(Node *left_,Node *right_,IntegerVector const &normal_, vector<int> const &conesInHyperplane_):
      normal(normal_),
      conesInHyperplane(conesInHyperplane_)
    {
      leftRight[0]=left_;
      leftRight[1]=right_;
    }
  };
  Node *root;
  bool hasBeenBuilt;
  static vector<SmallCone> buildPointers(vector<PolyhedralCone> const &v)
    {
    vector<SmallCone> ret;
    for(int i=0;i<v.size();i++)
      ret.push_back(SmallCone(&(v[i]),SymmetryGroup::identity(v[i].ambientDimension())));
    return ret;
    }
  void buildHyperplanes();
  Node *buildTree(PolyhedralCone const &region, vector<int> const &active);
  BSPTree(int n_, vector<SmallCone> const &theCones_, PolyhedralCone const *restrictingCone=0, bool doBuild=true);
  static bool isOnRightHandSide(IntegerVectorList const &omega, IntegerVector const &normal);
  static PolyhedralCone glue(PolyhedralCone A, PolyhedralCone B, IntegerVector const &normal);
  PolyhedralCone regionRek(PolyhedralCone const &C, IntegerVectorList const &omega, const Node *tree)const;
  PolyhedralCone region(IntegerVectorList const &omega)const;
  static void printRek(Node const *p);
  void print()const;
  int numberOfRegionsRek(Node const *v)const;
  int numberOfRegions()const;
  IntegerVector heuristic1(PolyhedralCone const &region, vector<int> const &active);
  IntegerVector heuristic2(PolyhedralCone const &region, vector<int> const &active);

  //Routines for avoiding BSP tree expansion:
  bool isInComplement(IntegerVector const &v)const;
  IntegerVector randomVectorInComplement()const;
  IntegerVector perturbationToVectorInComplement(IntegerVector const &u, IntegerVector const &normal)const;
  IntegerVector firstIntersectionNormal(IntegerVector const &u, IntegerVector const &p,vector<SmallCone>::const_iterator i)const;
  PolyhedralCone regionFast(IntegerVector const &v)const;
  /**
   * Computes the multiplicity at a generic perturbation of v in the hyperplane defined by normal.
   */
  int multiplicity(IntegerVector const &v, IntegerVector const &normal)const;
};

#endif /* BSPTREE_H_ */
