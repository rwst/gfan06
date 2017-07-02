#include "tropical2.h"

#include <stdlib.h>
#include <iostream>

#include "buchberger.h"
#include "division.h"
#include "tropical.h"
#include "wallideal.h"
#include "dimension.h"
#include "halfopencone.h"
#include "breadthfirstsearch.h"
#include "tropicalbasis.h"
#include "tropicalcurve.h"
#include "linalg.h"

#include "field_rationalfunctions2.h"

#include "timer.h"
#include "log.h"

static Timer tropicalPrincipalIntersectionTimer("Tropical principal intersection",1);

//////////////////////////////////////////////////////////////////////////
//Using stable intersections and rational functions

PolynomialSet buildH(PolynomialSet const &g, IntegerVector const &isGeneric, bool noHomog=false)
{
	int sum=isGeneric.sum();
	int n=isGeneric.size();
	PolynomialRing coefRing(g.getRing().getField(),sum);
	FieldRationalFunctions2 coefField(coefRing);
	PolynomialRing R(coefField,n-sum);

	PolynomialRing R3=R.withVariablesAppended("H");
	PolynomialSet ret(R);

	for(PolynomialSet::const_iterator p=g.begin();p!=g.end();p++)
	{
		Polynomial q(R);
		for(TermMap::const_iterator i=p->terms.begin();i!=p->terms.end();i++)
		{
			IntegerVector coefv(sum);int IC=0;
			IntegerVector expv(n-sum);int IE=0;
			for(int j=0;j<n;j++)
				if(isGeneric[j])
				{
					coefv[IC++]=i->first.exponent[j];
				}
				else
				{
					expv[IE++]=i->first.exponent[j];
				}
			q+=Term(coefField.polynomialToFraction(Polynomial(Term(i->second,Monomial(coefRing,coefv)))),Monomial(R,expv));
//			Polynomial p(coefRing);
//			p+=Term(i->m.c,Monomial(coefv));
		}
		ret.push_back(q);
	}
	if(noHomog)return ret;
	return ret.homogenization(R3);
	//We must make sure that the ideal remains homogeneous when turning some variables into parameters.
}

IntegerVector reconstruct(IntegerVector const &v1, IntegerVector const &isGeneric)
{
	int n=isGeneric.size();
	IntegerVector v2(n);
	int I=0;
	for(int i=0;i<n;i++)if(!isGeneric[i])v2[i]=v1[I++]-v1[v1.size()-1];
	return v2;
}

/**
 *
 */
IntegerVector nonTrivialTropismInner(PolynomialSet const &groebnerBasis)
{
	PolynomialRing R3=groebnerBasis.getRing().withVariablesAppended("H");
	PolynomialSet g=saturatedIdeal(groebnerBasis.homogenization(R3));
	int n=g.getRing().getNumberOfVariables()-1;
	int h=dimensionOfHomogeneitySpace(g)-1;
	int d=krullDimension(g)-1;

	//	greedy matroid algorithm
	IntegerVector isGeneric(n);
	int i=0;
	while(d>1)
	{
		assert(i<n);
		isGeneric[i]=1;

		bool OK;
		if(0)
		{	//Using rational functions as coefficients
			PolynomialSet r=buildH(groebnerBasis,isGeneric);
			OK=!containsMonomial(r);
		}
		else
		{	//Avoiding rational functions as coefficients
			OK=true;
			IntegerVectorList M;
			M.push_back(IntegerVector::allOnes(isGeneric.size())-isGeneric);
			M.push_back(IntegerVector::allOnes(isGeneric.size()));
			PolynomialSet G=groebnerBasis;
			buchberger(&G,MatrixTermOrder(M),true);
			for(PolynomialSet::const_iterator i=G.begin();i!=G.end();i++)
				if(!i->isZero())if(i->usedVariables().divides(isGeneric)){
					OK=false;break;}
		}

		if(!OK)
		{
			isGeneric[i]=0;
		}
		else
		{
			d--;
		}
		i++;
	}
//	debug<<"d="<<d<<"\n";
	if(d==1)
	{
	IntegerVectorList A=tropicalCurve(buildH(groebnerBasis,isGeneric,true),true);
		return reconstruct(A.front(),isGeneric);

		PolynomialSet r=buildH(groebnerBasis,isGeneric);
		buchberger(&r,LexicographicTermOrder());
		PolyhedralFan bergmanFan(r.numberOfVariablesInRing());
		PolynomialSet tropBasis=tropicalBasisOfCurve(r.numberOfVariablesInRing(),r,&bergmanFan,1);

		if(bergmanFan.dimensionOfLinealitySpace()==1)
		{
			IntegerVectorList rays=bergmanFan.getRays(2);
// At this point the lineality space of the ideal may have increased. A ray is sought outside the original lineality space.
			assert(!rays.empty());
			IntegerVector v1=*rays.begin();
			return reconstruct(v1,isGeneric);
		}
		else
		{
			assert(bergmanFan.dimensionOfLinealitySpace()==2);
			IntegerVectorList l=bergmanFan.conesBegin()->generatorsOfLinealitySpace();
			assert(l.size()==2);
			for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
			{
				if(!groebnerBasis.isHomogeneous(reconstruct(*i,isGeneric)))
					return reconstruct(*i,isGeneric);
			}
			pout<<l;
			pout<<groebnerBasis;
			assert(0);
		}
	}
	return IntegerVector(n);
}

IntegerVector nonTrivialTropism(PolynomialSet const &groebnerBasis)
{
	list<int> toKeep=groebnerBasis.multiDeHomogenizationToKeep();
	PolynomialSet ideal2=groebnerBasis.multiDeHomogenization();
	IntegerVector v=nonTrivialTropismInner(ideal2);
	IntegerVector ret(groebnerBasis.numberOfVariablesInRing());
	int I=0;
	for(list<int>::const_iterator i=toKeep.begin();i!=toKeep.end();i++,I++)
		ret[*i]=v[I];
	return ret;
}

static void startingConeError()
{
  fprintf(Stderr,"UNABLE TO COMPUTE STARTING CONE.\n");
  fprintf(Stderr,"THE STARTING CONE ALGORITHM IN GFAN IS BASED ON HEURISTICS WHICH HAVE FAILED ON THIS EXAMPLE.\n");
  assert(0);
}

PolynomialSet initialIdeal(PolynomialSet const &g, IntegerVector const &weight)
//Assume homogeneous
{
  PolynomialSet ret=g;
  WeightReverseLexicographicTermOrder T(weight);
  buchberger(&ret,T);
  return initialForms(ret,weight);
}

PolynomialSet initialIdealNonHomogeneous(PolynomialSet const &g, IntegerVector const &weight, bool allowSaturation)
{
	// We use Theorem 4 page 379 of [Cox, Little, O'Shea] (book 1) to compute the homogenization
	IntegerVector w(IntegerVector::allOnes(g.getRing().getNumberOfVariables()));
	WeightReverseLexicographicTermOrder T(w);
	PolynomialSet g2=g;
	buchberger(&g2,T,allowSaturation);
	PolynomialRing R2(g.getRing().getField(),g.getRing().getNumberOfVariables()+1);
	PolynomialSet g3=g2.homogenization(R2,&w);
	// We now compute the initial ideal using Proposition 5.2.3 in [Jensen] (PhD thesis).
	IntegerVector weight2=concatenation(weight,IntegerVector(1));
	WeightTermOrder/*WeightReverseLexicographicTermOrder*/ T2(weight2);// The proposition assumes that weight2 is tie-broken with a term order, but that is probably not necessary. However, the specifications of _this_ routine says that a GB must be returned.
	buchberger(&g3,T2,allowSaturation);
	PolynomialSet g4=initialFormsAssumeMarked(g3,weight2);
	return g4.deHomogenization();
}


Polynomial initialFormAssumeMarked(Polynomial const &p, IntegerVector const &weight)
{
  Polynomial r(p.getRing());
  IntegerVector markedExponent=p.getMarked().m.exponent;

  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    {
      IntegerVector dif=markedExponent-i->first.exponent;

      if(dot(dif,weight)==0)
	r+=Polynomial(Term(i->second,i->first));
    }
  r.mark(Monomial(p.getRing(),markedExponent));

  return r;
}


PolynomialSet initialFormsAssumeMarked(PolynomialSet const &groebnerBasis, IntegerVector const &weight)
{
  PolynomialRing theRing=groebnerBasis.getRing();
  PolynomialSet r(theRing);

  for(PolynomialSet::const_iterator i=groebnerBasis.begin();i!=groebnerBasis.end();i++)
    {
      r.push_back(initialFormAssumeMarked(*i,weight));
    }
  return r;
}


Polynomial initialForm(Polynomial const &p, IntegerVector const &weight)
{
  if(p.isZero())return p;
  int64 a=dotLong(p.terms.begin()->first.exponent,weight);
  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    {
      int64 b=dotLong(i->first.exponent,weight);
      if(b>a)a=b;
    }

  Polynomial r(p.getRing());
  bool ismarked=p.isMarked();
  IntegerVector markedExponent;
  if(ismarked)markedExponent=p.getMarked().m.exponent;

  bool markedFound=false;

  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    {
      if(dotLong(i->first.exponent,weight)==a)
	{
	  r+=Polynomial(Term(i->second,i->first));
	  if(ismarked)if((markedExponent-(i->first.exponent)).isZero())markedFound=true;
	}
    }
  if(markedFound)
    r.mark(Monomial(p.getRing(),markedExponent));

  return r;
}


PolynomialSet initialForms(PolynomialSet const &groebnerBasis, IntegerVector const &weight)
{
  PolynomialRing theRing=groebnerBasis.getRing();
  PolynomialSet r(theRing);
  if(theRing.getNumberOfVariables()!=weight.size())
    {
      cerr << "Error: Number of varaibles in polynomial ring "<<theRing.getNumberOfVariables()<< " length of weight vector " << weight.size() <<endl;
      assert(0);
    }

  for(PolynomialSet::const_iterator i=groebnerBasis.begin();i!=groebnerBasis.end();i++)
    {
      r.push_back(initialForm(*i,weight));
    }
  return r;
}



PolyhedralFan tropicalPrincipalIntersection(int n, PolynomialSet const &g, int linealitySpaceDimension)
{
  //return tropicalHyperSurfaceIntersection(n, g);////////////////////////////////////////
  log2 fprintf(Stderr,"Intersecting\n");
  log3 AsciiPrinter(Stderr).printPolynomialSet(g);

  TimerScope ts(&tropicalPrincipalIntersectionTimer);
  PolyhedralFan ret=PolyhedralFan::fullSpace(n);

  for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    {
      ret=refinement(ret,PolyhedralFan::bergmanOfPrincipalIdeal(*i),linealitySpaceDimension,true);
    }
  log2 fprintf(Stderr,"Done intersecting\n");
  return ret;
}



static PolynomialSet checkList(IntegerVectorList const &l, PolynomialSet const &groebnerBasis, PolynomialSet *fullNeighbourBasis, int h, bool &result, bool onlyCheckRays)
{
  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    {
      WeightReverseLexicographicTermOrder t(*i);
      log2 fprintf(Stderr,"Computing Gr\"obner basis with respect to:");
      log2 AsciiPrinter(Stderr).printVector(*i);
      log2 fprintf(Stderr,"\n");
      PolynomialSet h2=groebnerBasis;
      buchberger(&h2,t);
      log2 fprintf(Stderr,"Done computing Gr\"obner basis.\n");

      log3 AsciiPrinter(Stderr).printPolynomialSet(h2);
      PolynomialSet wall=initialFormsAssumeMarked(h2,*i);

      log3 AsciiPrinter(Stderr).printString("Initial ideal:\n");
      log3 AsciiPrinter(Stderr).printPolynomialSet(wall);

      int hdim2=dimensionOfHomogeneitySpace(wall);
      if(hdim2>h)
	{
	  if(!containsMonomial(wall))
	    {
	      log1 fprintf(Stderr,"Iterating recursively.\n");
	      //PolynomialSet initialIdeal=guessInitialIdealWithoutMonomial(wall,0);
	      PolynomialSet initialIdeal=guessInitialIdealWithoutMonomial(wall,fullNeighbourBasis,onlyCheckRays);

	      if(fullNeighbourBasis)
		{
		  //*fullNeighbourBasis=liftBasis(initialIdeal,h2);
		  *fullNeighbourBasis=liftBasis(*fullNeighbourBasis,h2);
		}


	      result=true;
	      return initialIdeal;
	    }
	}
    }
  result=false;
  return groebnerBasis;
}



IntegerVectorList perturbedRelativeInteriorTropism(PolynomialSet const &groebnerBasis)
{
	int n=groebnerBasis.getRing().getNumberOfVariables();
	// If T(I) is just the origin, we return the empty list
	PolynomialRing R(groebnerBasis.getRing().getField(),n+1);
	PolynomialSet g=groebnerBasis.homogenization(R);
	saturatedIdeal(g);
	if(krullDimension(g)==1)return IntegerVectorList();

	// If the homogeneity space is non-trivial find a basis of it and reduce dimension
	PolyhedralCone hom=homogeneitySpace(g);
	if(hom.dimension()>1)
	{
		FieldMatrix m=integerMatrixToFieldMatrix(rowsToIntegerMatrix(hom.getEquations(),n+1),Q);
		m.reduce();
		list<int> toKeep=m.pivotColumns();
		PolynomialRing R2=PolynomialRing(g.getRing().getField(),toKeep.size());

		PolynomialSet i2=g.embeddedInto(R2,&toKeep);

		IntegerVectorList rl=perturbedRelativeInteriorTropism(i2);

		IntegerVectorList eq;eq.push_back(IntegerVector::standardVector(n+1,n));
		PolyhedralCone C(IntegerVectorList(),eq,n+1);C=intersection(C,hom);C.canonicalize();
		IntegerVectorList ret=C.projection(n).generatorsOfLinealitySpace();
		for(IntegerVectorList::const_iterator j=rl.begin();j!=rl.end();j++)
		{
			IntegerVector temp(n+1);
			int I=0;
			for(list<int>::const_iterator i=toKeep.begin();i!=toKeep.end();i++,I++)
				temp[*i]=(*j)[I];
			ret.push_back(temp.subvector(0,n)-temp[n]*IntegerVector::allOnes(n));
		}
		return ret;
	}

	IntegerVector v=nonTrivialTropism(groebnerBasis);
	IntegerVectorList ret=perturbedRelativeInteriorTropism(initialIdealNonHomogeneous(groebnerBasis,v,true)); //At this point it would make sense to "dehomogenize w.r.t. v", so that v is not repeated in the recursion. We could also simply choose not to push v.
	ret.push_front(v);

	return ret;

/*	saturate
	if dim=0 return empty;
	compute homogeneity space
if homog>0
	reduce dimension of ring/ideal
	call recursively
	lift list
	return
	else
		return nontrivialtropism
	*/
}


PolynomialSet guessInitialIdealWithoutMonomial(PolynomialSet const &groebnerBasis, PolynomialSet *fullNeighbourBasis, bool onlyCheckRays) //ideal must be homogeneous
  // fullNeighbourBasis is set to a Groebner basis of the full ideal. The returned basis and fullNeighbourBasis have at least one termorder in common
{
	if(1)
	{
		IntegerVectorList w=perturbedRelativeInteriorTropism(groebnerBasis);
		MatrixTermOrder T(w);
		PolynomialSet g=groebnerBasis;
		buchberger(&g,T);//assuming that g is homogeneous
		PolynomialSet ig=g;
		for(IntegerVectorList::const_iterator i=w.begin();i!=w.end();i++)ig=initialForms(ig,*i);
		assert(fullNeighbourBasis);
		*fullNeighbourBasis=g;
		return ig;
	}
	//  log0 fprintf(Stderr,"A\n");
  assert(groebnerBasis.isValid());
  //  log0 fprintf(Stderr,"B\n");
  if(fullNeighbourBasis)
    {
      assert(fullNeighbourBasis->isValid());
    }
  //  log0 fprintf(Stderr,"C\n");

  int n=groebnerBasis.numberOfVariablesInRing();
  //  log0 fprintf(Stderr,"D\n");
  int h=dimensionOfHomogeneitySpace(groebnerBasis);
  //  log0 fprintf(Stderr,"E\n");
  int d=krullDimension(groebnerBasis);
  //  log0 fprintf(Stderr,"F\n");

  if(d==h)
    {
      if(fullNeighbourBasis)*fullNeighbourBasis=groebnerBasis;
      return groebnerBasis;
    }

#if 0
  // stable intersections/rational functions
  {
	  IntegerVectorList toCheck;
//	  pout<<groebnerBasis;


	  toCheck.push_back(nonTrivialTropism(groebnerBasis));
	  bool result;
	  pout<<toCheck<<n<<h<<d<<"\n";
	  PolynomialSet r=checkList(toCheck,groebnerBasis,fullNeighbourBasis,h,result, onlyCheckRays);
	  if(result)return r;
  }
  assert(0);
// TODO: When starting cone using stable intersections works, remove all dead code.
  #endif

  {
    //log2
    fprintf(Stderr,"Computing extreme rays.\n");
    //IntegerVectorList a;
    PolyhedralCone p=coneFromMarkedBasis(groebnerBasis);
    //PolyhedralCone p=PolyhedralCone(wallInequalities(groebnerBasis),a);
    IntegerVectorList extreme=p.extremeRays();
    log2 fprintf(Stderr,"Extreme rays of Groebner cone:\n");
    log2 AsciiPrinter(Stderr).printVectorList(extreme);

    bool result;
    PolynomialSet r=checkList(extreme,groebnerBasis,fullNeighbourBasis,h,result, onlyCheckRays);
    if(result)return r;
  }
  if(onlyCheckRays)startingConeError();
  PolyhedralFan f=PolyhedralFan::fullSpace(n);
  /*  for(int i=0;i<d-1;i++)
    {
      IntegerVector v(n);
      for(int j=0;j<n;j++)v[j]=rand()&1;
      IntegerVectorList a,b;
      b.push_back(v);
      PolyhedralFan F(n);
      F.insert(PolyhedralCone(a,b));
      f=refinement(f,F);
    }
  AsciiPrinter P(Stderr);
  f.print(&P);
  */
  int hypersurfacesToGo=groebnerBasis.size();
  for(PolynomialSet::const_iterator i=groebnerBasis.begin();i!=groebnerBasis.end();i++)
    {
      fprintf(Stderr,"Hypersurfaces to go:%i\n",hypersurfacesToGo--);
      fprintf(Stderr,"Max dimension: %i\n",f.getMaxDimension());
      debug<<"*i="<<*i<<"\n";
      f=refinement(f,PolyhedralFan::bergmanOfPrincipalIdeal(*i));
      f.removeAllExcept(3);

      IntegerVectorList l=f.getRelativeInteriorPoints();

      bool result;
      PolynomialSet r=checkList(l,groebnerBasis,fullNeighbourBasis,h,result, onlyCheckRays);
      if(result)return r;
    }
  startingConeError();
  return groebnerBasis;
}

static PolynomialSet checkListStably(IntegerVectorList const &l, PolynomialSet const &groebnerBasis, PolynomialSet *fullNeighbourBasis, int h, bool &result, bool onlyCheckRays)
{
  debug<< "Checklist called on"<<groebnerBasis;
  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    {
      WeightReverseLexicographicTermOrder t(*i);
      log2 fprintf(Stderr,"Taking initial forms with respect to:");
      log2 AsciiPrinter(Stderr).printVector(*i);
      log2 fprintf(Stderr,"\n");
      PolynomialSet h2=groebnerBasis;
      log2 fprintf(Stderr,"Done computing Gr\"obner basis.\n");

      log3 AsciiPrinter(Stderr).printPolynomialSet(h2);
      PolynomialSet wall=initialForms(h2,*i);

      log3 AsciiPrinter(Stderr).printString("Initial ideal:\n");
      log3 AsciiPrinter(Stderr).printPolynomialSet(wall);

      int hdim2=dimensionOfHomogeneitySpace(wall);
      if(hdim2>h)
        {
          if(nonEmptyStableIntersection(wall))
            {
              log1 fprintf(Stderr,"Iterating recursively.\n");
              //PolynomialSet initialIdeal=guessInitialIdealWithoutMonomial(wall,0);
              PolynomialSet initialIdeal=guessInitialIdealWithoutMonomialStably(wall,fullNeighbourBasis,onlyCheckRays);

              if(fullNeighbourBasis)
                {
                  //*fullNeighbourBasis=liftBasis(initialIdeal,h2);
//                  *fullNeighbourBasis=liftBasis(*fullNeighbourBasis,h2);
                  *fullNeighbourBasis=groebnerBasis;
                  fullNeighbourBasis->copyMarkings(initialIdeal);
                }


              result=true;
              return initialIdeal;
            }
        }
    }
  result=false;
  return groebnerBasis;
}

PolynomialSet guessInitialIdealWithoutMonomialStably(PolynomialSet const &groebnerBasis, PolynomialSet *fullNeighbourBasis, bool onlyCheckRays) //ideal must be homogeneous
  // fullNeighbourBasis is set to a Groebner basis of the full ideal. The returned basis and fullNeighbourBasis have at least one termorder in common
{
  int n=groebnerBasis.numberOfVariablesInRing();
  int h=dimensionOfHomogeneitySpace(groebnerBasis);
  int d=n-groebnerBasis.size();//krullDimension(groebnerBasis);

  debug<</*"d"<<d<<*/"h"<<h<<"n"<<n<<"\n";


  if(d==h)
    {
      if(fullNeighbourBasis)*fullNeighbourBasis=groebnerBasis;
      return groebnerBasis;
    }

  {
    log2 fprintf(Stderr,"Computing extreme rays.\n");
    //IntegerVectorList a;
    PolyhedralCone p=coneFromMarkedBasis(groebnerBasis);
    //PolyhedralCone p=PolyhedralCone(wallInequalities(groebnerBasis),a);
    IntegerVectorList extreme=p.extremeRays();
    log2 fprintf(Stderr,"Extreme rays of Groebner cone:\n");
    log2 AsciiPrinter(Stderr).printVectorList(extreme);

    bool result;
    PolynomialSet r=checkListStably(extreme,groebnerBasis,fullNeighbourBasis,h,result, onlyCheckRays);
    if(result)return r;
  }
  if(onlyCheckRays)startingConeError();

  PolyhedralFan f=PolyhedralFan::fullSpace(n);

  int hypersurfacesToGo=groebnerBasis.size();
  for(PolynomialSet::const_iterator i=groebnerBasis.begin();i!=groebnerBasis.end();i++)
    {
      fprintf(Stderr,"Hypersurfaces to go:%i\n",hypersurfacesToGo--);
      fprintf(Stderr,"Max dimension: %i\n",f.getMaxDimension());
      f=refinement(f,PolyhedralFan::bergmanOfPrincipalIdeal(*i));
      f.removeAllExcept(3);

      IntegerVectorList l=f.getRelativeInteriorPoints();

      bool result;
      PolynomialSet r=checkListStably(l,groebnerBasis,fullNeighbourBasis,h,result, onlyCheckRays);
      if(result)return r;
    }
  startingConeError();
  return groebnerBasis;
}
