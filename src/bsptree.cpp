/*
 * bsptree.cpp
 *
 *  Created on: Feb 5, 2011
 *      Author: anders
 */

#include "linalg.h"
#include "bsptree.h"
#include "lp.h"
#include "matrix.h"
#include "log.h"

  void BSPTree::buildHyperplanes()
  {
    set<IntegerVector> temp;
    IntegerVector v;
    for(vector<SmallCone>::const_iterator i=theCones.begin();i!=theCones.end();i++)
      {
        i->assignEquation(v);
        temp.insert(v);
      }
    hyperplanes=vector<IntegerVector>(temp.size());
    int I=0;
    for(set<IntegerVector>::const_iterator i=temp.begin();i!=temp.end();i++,I++)
      hyperplanes[I]=*i;
  }

  /* Figure out what the right assumptions are. These seem not to be right:
   * A is assumed to be full-dimensional.
   * B is assumed to be of codimension 1 with known facets and implied equations.
   */
  static bool isIntersectionOfCodimensionOne(PolyhedralCone const &A, BSPTree::SmallCone const &B)
  {

    IntegerVectorList strictInequalities=A.getHalfSpaces();
    IntegerVectorList strictInequalities2=B.getHalfSpaces();
    for(IntegerVectorList::const_iterator i=strictInequalities2.begin();i!=strictInequalities2.end();i++)strictInequalities.push_back(*i);
    if(B.getEquations().size()!=1)
      {
assert(0);        //debug<<B;
      }
    assert(B.getEquations().size()==1);
    strictInequalities.push_front(B.getEquations().front());
    IntegerVector ei(strictInequalities.size());ei[0]=1;
//    bool originalReturnValue=hasInteriorPoint(strictInequalities, false, &ei);


    //////////////////
    IntegerVectorList inequalities;
      IntegerVectorList::const_iterator i=strictInequalities.begin();i++;
      for(;i!=strictInequalities.end();i++)
        inequalities.push_back(concatenation(-IntegerVector::standardVector(1,0),*i));
      IntegerVectorList equations;
      equations.push_back(concatenation(IntegerVector(1),strictInequalities.front()));
      bool fastReturnValue=hasHomogeneousSolution(A.ambientDimension()+1,inequalities,equations);

  //    assert(originalReturnValue==fastReturnValue);
      /////////////////
    return fastReturnValue;

//    return originalReturnValue;
  }
  static bool isInClosedHalfSpace(BSPTree::SmallCone const &A, IntegerVector const &normal)
  {
    //WE ALSO NEED TO CHECK LINEALITYSPACE!
    IntegerVectorList temp2=A.generatorsOfLinealitySpace();
    for(IntegerVectorList::const_iterator i=temp2.begin();i!=temp2.end();i++)
      if(dotLong(normal,*i)!=0)return false;



    IntegerVectorList temp=A.extremeRays();
    for(IntegerVectorList::const_iterator i=temp.begin();i!=temp.end();i++)
      if(dotLong(normal,*i)<0)return false;
    return true;
  }


  class ProgressPrinter
  {
public:
    static int depth;
    static char s[16];
    static int timeToPrint;
  public:
    ProgressPrinter()
    {
      if(depth>0 && depth<17)s[depth-1]++;
      depth++;
      if(depth<17)s[depth-1]=0;

      timeToPrint++;
      if(depth<16)
//      if(!(timeToPrint&1))
        {
          log1 debug<<"Progress:";
          log1 for(int i=0;(i<depth-1)&&(i<16);i++)
              debug<<((s[i]-1)?"1":"0");
          log1 debug<<"\n";
        }
    }
    ~ProgressPrinter()
    {
      depth--;
    }
  };
  int ProgressPrinter::depth;
  char ProgressPrinter::s[16];
  int ProgressPrinter::timeToPrint;

  IntegerVector BSPTree::heuristic1(PolyhedralCone const &region, vector<int> const &active)
    {
    int index=active[rand()%active.size()];
    assert(theCones[index].getEquations().size()==1);
    IntegerVector normal=theCones[index].getEquations().front();

    {
      map<IntegerVector,int> M;
      for(vector<int>::const_iterator i=active.begin();i!=active.end();i++)
        {
          M[theCones[*i].getEquations().front()]++;
        }
      int highScore=0;
      IntegerVector const *p=0;
      for(map<IntegerVector,int>::const_iterator i=M.begin();i!=M.end();i++)
        {
          if(highScore<i->second)
            {
              highScore=i->second;
              p=&(i->first);
            }
        }
      normal=*p;
    }
    return normal;
    }

  IntegerVector BSPTree::heuristic2(PolyhedralCone const &region, vector<int> const &active)
  {
    if(active.size()<1000)return heuristic1(region,active);

    int bestScore=1000;
    IntegerVector bestNormal=theCones[active[0]].getEquations().front();

    for(int i=0;i<100;i++)
      {
        int index=active[rand()%active.size()];
        int inPlane=0;
        int right=0;
        int left=0;

        IntegerVector normal=theCones[index].getEquations().front();

        IntegerVectorList temp;temp.push_back(normal);
        IntegerVectorList empty;
        PolyhedralCone rightRegion=intersection(region,PolyhedralCone(temp,empty));
        PolyhedralCone leftRegion=intersection(region,PolyhedralCone(temp,empty).negated());
        PolyhedralCone hyperPlane=intersection(region,PolyhedralCone(empty,temp));


        for(int j=0;j<100;j++)
          {
//          int r=
              int rindex=active[rand()%active.size()];
          if(dependent(normal,theCones[rindex].getEquations().front()))inPlane++;
          else
            {
            if(isIntersectionOfCodimensionOne(leftRegion,theCones[rindex]))left++;
            if(isIntersectionOfCodimensionOne(rightRegion,theCones[rindex]))right++;
            }
          }
        int score=left+right;

        debug<<inPlane<<" "<<left<<" "<<right<<" "<<left+right<<"\n";

        if(score<bestScore)
          {
            bestScore=score;
            bestNormal=normal;
          }
      }
      return bestNormal;
//    return heuristic1(region,active);
  }


  BSPTree::Node *BSPTree::buildTree(PolyhedralCone const &region, vector<int> const &active)
  {
    ProgressPrinter PP;
  //  depth++;
  //  debug<<"Size:"<<(int)active.size()<<"\n";
    if(active.size()==0)return 0;

    IntegerVector normal=heuristic1(region,active);

    vector<int> conesInHyperplane;
    vector<int> conesLeft;
    vector<int> conesRight;
    IntegerVectorList temp;temp.push_back(normal);
    IntegerVectorList empty;
    PolyhedralCone rightRegion=intersection(region,PolyhedralCone(temp,empty));
    PolyhedralCone leftRegion=intersection(region,PolyhedralCone(temp,empty).negated());
    PolyhedralCone hyperPlane=intersection(region,PolyhedralCone(empty,temp));
    int counter=0;
    for(vector<int>::const_iterator i=active.begin();i!=active.end();i++)
      {
        counter++;
        log1 if(!(counter&(4096*4-1)))debug<<counter<<"\n";
        /*
         * Cases:
         * The intersection of the cone with the relative interior of the left region is of codimension 1
         * The intersection of the cone with the relative interior of the right region is of codimension 1
         * None of these are the case, then the cone is contained in the chosen hyperplane, and should be added to the node.
         */

//        PolyhedralCone iHyperplane=intersection(theCones[*i],hyperPlane);
//        iHyperplane.findImpliedEquations();
//        if(iHyperplane.dimension()==n-1)
        if(dependent(normal,theCones[*i].getEquations().front()))
          conesInHyperplane.push_back(*i);
        else
          {
            bool A=isInClosedHalfSpace(theCones[*i],-normal);
            bool B=isInClosedHalfSpace(theCones[*i],normal);
//            if(isInClosedHalfSpace(theCones[*i],-normal)||!isInClosedHalfSpace(theCones[*i],normal))
            if(A||!B)
            if(isIntersectionOfCodimensionOne(leftRegion,theCones[*i]))
              conesLeft.push_back(*i);
//            if(isInClosedHalfSpace(theCones[*i],normal)||!isInClosedHalfSpace(theCones[*i],-normal))
            if(B||!A)
            if(isIntersectionOfCodimensionOne(rightRegion,theCones[*i]))
              conesRight.push_back(*i);
          }
      }
    if(PP.depth<10+15)debug<<PP.depth<<":"<<(int)conesInHyperplane.size()<<" "<<(int)conesLeft.size()<<" "<<(int)conesRight.size()<<"\n";
    return new Node(
        buildTree(leftRegion,conesLeft),
        buildTree(rightRegion,conesRight),
        normal,conesInHyperplane);
  }
  BSPTree::BSPTree(int n_, vector<BSPTree::SmallCone> const &theCones_, PolyhedralCone const *restrictingCone, bool doBuild):
    theCones(theCones_),
    n(n_),
    hasBeenBuilt(false)
    {
      buildHyperplanes();
      vector<int> active(theCones.size());
      for(int i=0;i<active.size();i++)active[i]=i;
      if(doBuild)
        {
          if(restrictingCone)
            root=buildTree(*restrictingCone,active);
          else
            root=buildTree(PolyhedralCone(n),active);
          hasBeenBuilt=true;
        }
      else
        root=0;
    }
  bool BSPTree::isOnRightHandSide(IntegerVectorList const &omega, IntegerVector const &normal)
  {
    for(IntegerVectorList::const_iterator i=omega.begin();i!=omega.end();i++)
      {
        if(dotLong(*i,normal)>0)return true;
        if(dotLong(*i,normal)<0)return false;
      }
    for(int i=0;i<normal.size();i++)
      {
        if(normal[i]>0)return true;
        if(normal[i]<0)return false;
      }
    assert(0);
    return true;
  }
  PolyhedralCone BSPTree::glue(PolyhedralCone A, PolyhedralCone B, IntegerVector const &normal)
  {
    A.findFacets();
    B.findFacets();
    IntegerVectorList inequalitiesA=A.getHalfSpaces();
    IntegerVectorList inequalitiesB=B.getHalfSpaces();
    set<IntegerVector> inequalities;
    for(IntegerVectorList::const_iterator i=inequalitiesA.begin();i!=inequalitiesA.end();i++)
      if(!dependent(*i,normal))inequalities.insert(normalized(*i));
    for(IntegerVectorList::const_iterator i=inequalitiesB.begin();i!=inequalitiesB.end();i++)
      if(!dependent(*i,normal))inequalities.insert(normalized(*i));
    IntegerVectorList empty;
    IntegerVectorList inequalities2;for(set<IntegerVector>::const_iterator i=inequalities.begin();i!=inequalities.end();i++)inequalities2.push_back(*i);
    PolyhedralCone ret(inequalities2,empty,normal.size(),PCP_facetsKnown|PCP_impliedEquationsKnown);
    return ret;
//    ret.canonicalize();

#if 0
    //Assumes that A and B are full dimensional.
    A.canonicalize();
    B.canonicalize();
/*    debug<<"Gluing:\n";
    debug<<normal<<"\n";
    debug<<A;
    debug<<B;
*/
  /*  {
      PolyhedralCone in=intersection(A,B);
      IntegerVector v=in.getRelativeInteriorPoint();
      PolyhedralCone FA=A.faceContaining(v);
      PolyhedralCone FB=B.faceContaining(v);

      FA.canonicalize();
      FB.canonicalize();
      debug<<"FA"<<FA;
      debug<<"FB"<<FB;
      assert(!(FA<FB));
      assert(!(FB<FA));
    }
*/


    IntegerVectorList inequalitiesA=A.getHalfSpaces();
    IntegerVectorList inequalitiesB=B.getHalfSpaces();
    IntegerVectorList inequalities;
    for(IntegerVectorList::const_iterator i=inequalitiesA.begin();i!=inequalitiesA.end();i++)
      if(!dependent(*i,normal))inequalities.push_back(*i);
    for(IntegerVectorList::const_iterator i=inequalitiesB.begin();i!=inequalitiesB.end();i++)
      if(!dependent(*i,normal))inequalities.push_back(*i);

    //debug<<inequalities;
    IntegerVectorList empty;
    PolyhedralCone ret(inequalities,empty,normal.size());
    ret.canonicalize();
    //debug<<"RET"<<ret;
/*
    {
      IntegerVectorList generatorsA=A.extremeRays();
      IntegerVectorList generatorsB=B.extremeRays();
      IntegerVectorList generatorsAL=A.generatorsOfLinealitySpace();
      IntegerVectorList generatorsBL=B.generatorsOfLinealitySpace();
      for(IntegerVectorList::const_iterator i=generatorsB.begin();i!=generatorsB.end();i++)generatorsA.push_back(*i);
      for(IntegerVectorList::const_iterator i=generatorsBL.begin();i!=generatorsBL.end();i++)generatorsAL.push_back(*i);
      PolyhedralCone ret2=PolyhedralCone::givenByRays(generatorsA,generatorsAL,A.ambientDimension());
      ret2.canonicalize();
      debug<<"RET2"<<ret2;
      assert(ret.contains(ret2));
      assert(ret2.contains(ret));
      assert(!(ret2<ret));
      assert(!(ret<ret2));
//      return ret2;
//      assert(ret==ret2);
    }*/
    return ret;
#endif
    }

  PolyhedralCone BSPTree::regionRek(PolyhedralCone const &C, IntegerVectorList const &omega, const Node *tree)const
  {
    if(tree==0)return C;
    IntegerVectorList empty;
    IntegerVectorList temp;temp.push_back(tree->normal);
    PolyhedralCone hyperPlane(empty,temp);
    PolyhedralCone rightHalfSpace(temp,empty);//HERE
    PolyhedralCone leftHalfSpace=rightHalfSpace.negated();
    int sideWeAreOn=isOnRightHandSide(omega,tree->normal);
    IntegerVectorList tempIneq=C.getHalfSpaces();tempIneq.push_back(sideWeAreOn?tree->normal:-tree->normal);
    PolyhedralCone tempIntersection(tempIneq,empty,n,PCP_impliedEquationsKnown);
//    PolyhedralCone A=regionRek(intersection(C,sideWeAreOn?rightHalfSpace:leftHalfSpace),omega,tree->leftRight[sideWeAreOn]);
    PolyhedralCone A=regionRek(tempIntersection,omega,tree->leftRight[sideWeAreOn]);

    bool isANormal=false;
    A.findFacets();
    IntegerVectorList normals=A.getHalfSpaces();
    for(IntegerVectorList::const_iterator i=normals.begin();i!=normals.end();i++)
      if(dependent(tree->normal,*i)){isANormal=true;break;}
    if(!isANormal)return A;


    PolyhedralCone iA(normals,temp,n,PCP_impliedEquationsKnown);
//    PolyhedralCone iA=intersection(A,hyperPlane);
//    iA.findImpliedEquations();


//    if(iA.dimension()==n-1)
      {
      IntegerVector v=iA.getRelativeInteriorPoint();
      bool isContained=false;
      for(vector<int>::const_iterator i=tree->conesInHyperplane.begin();i!=tree->conesInHyperplane.end();i++)
        if(theCones[*i].contains(v)){isContained=true;break;}
      if(!isContained)
        {
          IntegerVectorList omega2;
          omega2.push_back(v);omega2.push_back(sideWeAreOn?-(tree->normal):tree->normal);//HERE
          IntegerVectorList tempIneq2=C.getHalfSpaces();tempIneq2.push_back(sideWeAreOn?-tree->normal:tree->normal);
          PolyhedralCone tempIntersection2(tempIneq2,empty,n,PCP_impliedEquationsKnown);
          PolyhedralCone B=regionRek(tempIntersection2,omega2,tree->leftRight[1-sideWeAreOn]);
//          PolyhedralCone B=regionRek(intersection(C,sideWeAreOn?leftHalfSpace:rightHalfSpace),omega2,tree->leftRight[1-sideWeAreOn]);
          return glue(A,B,tree->normal);
        }
      }
    return A;
  }
  PolyhedralCone BSPTree::region(IntegerVectorList const &omega)const
  {
    assert(hasBeenBuilt);
    PolyhedralCone C=PolyhedralCone(n);
    return regionRek(C,omega,root);
  }
  void BSPTree::printRek(Node const *p)
  {
    if(p)
      {
        debug<<p->normal;
        debug<<"[";
        printRek(p->leftRight[1]);
        debug<<":";
        printRek(p->leftRight[0]);
        debug<<"]\n";
      }

  }
  void BSPTree::print()const
  {
    printRek(root);
  }
  int BSPTree::numberOfRegionsRek(Node const *v)const
  {
    if(v==0)return 1;
    return numberOfRegionsRek(v->leftRight[0])+numberOfRegionsRek(v->leftRight[1]);
  }
  int BSPTree::numberOfRegions()const
  {
    return numberOfRegionsRek(root);
  }

  static bool isPerturbedDotProductPositive(IntegerVector const &p, IntegerVector const &ineq)
    {
      if(dotLong(p,ineq)>0)return true;
      if(dotLong(p,ineq)<0)return false;
      for(int i=0;i<p.size();i++)
        {
          if(ineq[i]>0)return true;
          if(ineq[i]<0)return false;
        }
      assert(ineq.isZero());
      return false;
    }

  /**
   * Given vectors u,p and hyperplanes h1 h2, decide if the ray going from u to a lex perturbed p intersects h1 before h2.
   */
  static bool before(IntegerVector const &u, IntegerVector const &p, IntegerVector const &h1, IntegerVector const &h2)
  {
    /*
     * The ray is (after scaling) parameterized as follows
     * v(t)=tu+b+eps*e_1+...+eps^n*e_n
     * where t goes from infinity through 0 to minus infinity.
     *
     * The times for intersection are
     * t_i=(hi*b+eps*hi_1+...+eps^n*hi_n)/(-u*hi)
     * WLOG we assume that u*hi>0
     *
     * Hence we check if
     * (h1*b+eps*h1_1+...+eps^n*h1_n)*(u*h2) < (h2*b+eps*h2_1+...+eps^n*h2_n)*(u*h1)
     */
/*    int64 uh1=dotLong(u,h1);
    int64 uh2=dotLong(u,h2);
    int64 ph1=dotLong(p,h1);
    int64 ph2=dotLong(p,h2);
*/
    int64 uh1,uh2,ph1,ph2;
    dotLong4(u,p,h1,h2,uh1,uh2,ph1,ph2);
    bool flag=(uh1<0)^(uh2<0);
/*    if(uh1<0)
      {
        uh1=-uh1;
        h1=-h1;
      }
    if(uh2<0)
      {
        uh2=-uh2;
        h2=-h2;
      }
    if(dotLong(h1,p)*uh2>dotLong(h2,p)*uh1)return false^flag;
    if(dotLong(h1,p)*uh2<dotLong(h2,p)*uh1)return true^flag;
*/
    if(ph1*uh2>ph2*uh1)return false^flag;
    if(ph1*uh2<ph2*uh1)return true^flag;
    for(int i=0;i<u.size();i++)
      {
        if(h1[i]*uh2>h2[i]*uh1)return false^flag;
        if(h1[i]*uh2<h2[i]*uh1)return true^flag;
      }
    if(flag)assert((h1+h2).isZero());
    else assert((h1-h2).isZero());
    return false;
  }

/**
 * Checks if the given vector is in the complement of the stored cones.
 */
  bool BSPTree::isInComplement(IntegerVector const &v)const
  {
    for(vector<SmallCone>::const_iterator i=theCones.begin();i!=theCones.end();i++)
      if(i->contains(v))return false;
    return true;
  }

  /**
   * Returns random vector in complement.
   */
  IntegerVector BSPTree::randomVectorInComplement()const
  {
    while(1)
      {
        IntegerVector ret=randomIntegerVector(n,100);

        if(isInComplement(ret))return ret;
      }
    return IntegerVector();
  }


  /**
   * Given a vector u and a vector n, such that none of the stored cones contain u+epsilon n for epsilon>0 sufficient small,
   * find a vector ret in the span of u and n such that the open line segment between u and ret does intersect any of the
   * stored cones.
   */
  IntegerVector BSPTree::perturbationToVectorInComplement(IntegerVector const &u, IntegerVector const &normal)const
  {
/*    if(theCones.size()==0)return u+normals;
    vector<SmallCone>::const_iterator best=theCones.end();
    for(vector<SmallCone>::const_iterator i=theCones.begin();i!=theCones.end();i++)
      {
        IntegerVector equation=i->getEquation().front();
        if(dotLong(equation,u)==0)continue;

      }

      */
    //Fix this later

    bool found=false;
    IntegerVector const *current=0;
    for(vector<IntegerVector>::const_iterator i=hyperplanes.begin();i!=hyperplanes.end();i++)
      {
        if(dotLong(*i,u)){
            current=&*i;
            found=true;
        }
      }
    if(!found)
      {
        return normal;
      }
    for(vector<IntegerVector>::const_iterator i=hyperplanes.begin();i!=hyperplanes.end();i++)
      if(dotLong(*i,u))
        if(before(u, normal, *i, *current))current=&*i;

    int s=1;
    int64 ud=dotLong(*current,u);
    IntegerVector ret=u+normal;
    for(int i=0;i<30;i++)
      {
        if(ud<0 && dotLong(*current,ret)<0)return ret;
        if(ud>0 && dotLong(*current,ret)>0)return ret;
        s*=2;
        ret+=s*u;
      }
    assert(0);
    return 100*u+normal;
  }


  /**
   * Given a subset of the stored cones, a vector u in the complement and a point p, find the first codimension one intersecting
   * the ray starting at u and passing through a lex perturbation of p.
   */
  IntegerVector BSPTree::firstIntersectionNormal(IntegerVector const &u, IntegerVector const &p, vector<SmallCone>::const_iterator i)const
  {
    IntegerVector ret;
    bool found=false;
    IntegerVector equation;
//    for(vector<SmallCone>::const_iterator i=theCones.begin();i!=theCones.end();i++)
      for(;i!=theCones.end();i++)
      {

        i->assignEquation(equation);
      //  IntegerVectorList equations=i->getEquations();
      //  assert(equations.size()==1);
      //  IntegerVector equation=equations.front();
        // check if equation is better that current
        if(dotLong(equation,u))
        if(!found || before(u,p,equation,ret))
          {
            bool doIntersect=true;
            int doint=i->doIntersectNonPerturbed(u,p);
            if(doint==2)goto noTestNeeded;
            if(doint==0)continue;
            //check if ray intersects cone
            {
                IntegerVectorList inequalities=i->getHalfSpaces();
            for(IntegerVectorList::const_iterator j=inequalities.begin();j!=inequalities.end();j++)
              {
                int64 ju=dotLong(*j,u);


                if( (ju==0)&& !isPerturbedDotProductPositive(p,*j)||
                    (ju>0)&&before(u,p,*j,equation)||
                    (ju<0)&&before(u,p,equation,*j))
                  {
                    doIntersect=false;
                    break;
                  }
              }
            }
            noTestNeeded:
            if(doIntersect)
              {
/*                if(ret.size())
                for(int k=0;k<12;k++)
                  {
                    IntegerVector temp=k*p+(10-k)*u;
                    D((int)dotLong(ret,temp));
                    D((int)dotLong(equation,temp));
                  }*/
                ret=equation;
                found=true;
            }
          }
      }
    assert(found);
    return ret;
  }


  /**
   * Given a vector u not contained in any of the stored cones find the closure of the connected component of the complement
   * of the stored cones containing u.
   */
  PolyhedralCone BSPTree::regionFast(IntegerVector const &u)const
  {
if(0)    {
      static int abort;
      if(!abort)D(abort);
      abort++;
      if(abort>1000)assert(0);
    }
    assert(u.size()==n);
    PolyhedralCone ret=(n);
    IntegerVectorList ineq;
    IntegerVectorList ineqNeg;
    for(vector<SmallCone>::const_iterator i=theCones.begin();i!=theCones.end();)
      {
        bool skip=false;
#if 0
        IntegerVectorList gen;i->appendGeneratorsOfLinealitySpace(gen);
        IntegerVectorList rays;i->appendExtremeRays(rays);
//        IntegerVectorList ineq=ret.getHalfSpaces();
        for(IntegerVectorList::const_iterator j=ineq.begin();j!=ineq.end();j++)
          {
            bool isCert=true;
            for(IntegerVectorList::const_iterator i=gen.begin();i!=gen.end();i++)
              if(dotLong(*i,*j)!=0){isCert=false;goto leave;}
            for(IntegerVectorList::const_iterator i=rays.begin();i!=rays.end();i++)
              if(dotLong(*i,*j)>0){isCert=false;goto leave;}
        leave:
          if(isCert){skip=true;break;}
          }
#else
        for(IntegerVectorList::const_iterator j=ineqNeg.begin();j!=ineqNeg.end();j++)
          {
            if(i->doesSatisfyInequalityExpensive(*j)){skip=true;break;}
          }
#endif
//          D(skip);

          if(!skip && isIntersectionOfCodimensionOne(ret,*i))
          {
            PolyhedralCone B=intersection(ret,PolyhedralCone(i->getHalfSpaces(),i->getEquations(),n));
            IntegerVector p=B.getRelativeInteriorPoint();
            IntegerVector w=firstIntersectionNormal(u,p,i);
            if(dotLong(w,u)<0)w=-w;
            IntegerVectorList temp;temp.push_back(w);
            ret=intersection(ret,PolyhedralCone(temp,IntegerVectorList(),n));
            ineq.push_back(w);
            ineqNeg.push_back(-w);
            assert(ret.contains(u));
          }
          else
            i++;
      }
    return ret;
  }

  int BSPTree::multiplicity(IntegerVector const &v, IntegerVector const &normal)const
  {
	  int ret=0;
	  IntegerVectorList temp;temp.push_back(normal);
	  IntegerMatrix t=rowsToIntegerMatrix(temp);
	  FieldMatrix M=integerMatrixToFieldMatrix(t,Q);
	  IntegerVectorList l=fieldMatrixToIntegerMatrixPrimitive(M.reduceAndComputeKernel()).getRows();
	  l.push_front(v);
	  for(vector<SmallCone>::const_iterator i=theCones.begin();i!=theCones.end();i++)
		  if(i->contains(v))
		  {
			  ret+=i->multiplicity();
		  }

	  return ret;
  }
