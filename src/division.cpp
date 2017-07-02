#include "division.h"
#include "printer.h"
#include "wallideal.h"
#include "lp.h"

#include "timer.h"

static Timer divisionTimer("Division",1);
static Timer divisionTimer1("Division1",1);
static Timer divisionTimer2("Division2",1);
static Timer divisionTimer3("Division3",1);
static Timer divisionTimer4("Division4",1);


Polynomial division1(Polynomial p, PolynomialSet const &l, TermOrder const &termOrder, PolynomialSet *q, bool earlyExit=false);



typedef map<Monomial,Polynomial,TermMapCompare> ReductionCache;



Polynomial division(Polynomial p, PolynomialSet const &l, TermOrder const &termOrder, PolynomialSet *q, bool earlyExit)
{
  return division1(p,l,termOrder,q,earlyExit);
}

Polynomial smartDivision(Polynomial p, PolynomialSet l, TermOrder const &termOrder)
{
  Polynomial r(p.getRing());
  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    {
      r+=division(Term(i->second,Monomial(p.getRing(),i->first.exponent)),l,termOrder);
    }
  return r;
}


IntegerVector termorderWeight(PolynomialSet const &g)
{
  IntegerVector ret=relativeInteriorPoint(g.getRing().getNumberOfVariables(),fastNormals(wallInequalities(g)),0);

  //  fprintf(stderr,"WEIGHT");
  //AsciiPrinter(Stderr).printVector(ret);
  //fprintf(stderr,"\n");
  return ret;
}

Polynomial division1(Polynomial p, PolynomialSet const &l, TermOrder const &termOrder, PolynomialSet *q, bool earlyExit)
{
  PolynomialRing theRing=p.getRing();
  TimerScope ts(&divisionTimer);

  //  WeightReverseLexicographicTermOrder termOrder2(termorderWeight(l));//REMOVE ME ?? JAN 2009

//{static int i;i++;if(i>2600000)pout<<p<<"\n";}
  if(q)
    {
      *q=PolynomialSet(theRing);
      for(PolynomialSet::const_iterator i=l.begin();i!=l.end();i++)
	q->push_back(Polynomial(p.getRing()));
    }
  Polynomial r(p.getRing());

//  pout<<l;
//  pout<<"\n\n";
  while(!p.isZero())
    {
//	  pout<<"p:"<<p<<"\n";
      //      AsciiPrinter(Stderr).printPolynomial(p);
      //     fprintf(Stderr,"Number Of terms: %i\n",p.terms.size());
      p.mark(termOrder);

      Term initial=p.getMarked();

      PolynomialSet::const_iterator i;
      PolynomialSet::const_iterator iBest=l.end();
      int bestLength=-1;

      PolynomialSet::iterator j;
      PolynomialSet::iterator jBest;
      if(q){j=q->begin();jBest=j;}

{
    	  {
    		  static int k;
    		  k++;
//    		  if((k%1000)==0)pout<<"Kk"<<k<<"\n";
    	  }
	TimerScope ts(&divisionTimer2);
	  for(i=l.begin();i!=l.end();i++)
	    {
	      if(i->getMarked().m.exponent.divides(initial.m.exponent)){break;pout<<"..."<<*i<<"\n";/*break;*/int length=i->numberOfTerms();if((bestLength==-1)||(length<bestLength)){bestLength=length;iBest=i;if(q)jBest=j;}}//break;
	      if(q)j++;
	    }
      }
      if(bestLength==0 &&!l.empty())
      {
    	  //fprintf(stderr,"%i\n",bestLength);
    	  pout<<l;
    	  assert(0);
      }

//      i=iBest;j=jBest;
      {
	TimerScope ts(&divisionTimer3);
	if(i!=l.end())
	  {
	    Term s(-initial.c*i->getMarked().c.inverse(),Monomial(p.getRing(),initial.m.exponent-i->getMarked().m.exponent));
//	    pout<<"OLD "<<p.numberOfTerms()<<" ";
	    p.madd(s,*i);

//	    pout<<"NEW "<<p.numberOfTerms()<<" I "<<i->numberOfTerms()<<"\n";

	    if(q)*j+=Polynomial(s);
	  }
	else
	  {
	    TimerScope ts(&divisionTimer4);
	    p-=initial;
	    r+=initial;
	    if(earlyExit)return r+p;
	  }
      }
    }
  return r;
}

#if 0
Polynomial divisionLift(Polynomial p, PolynomialSet l, PolynomialSet lLift, TermOrder const &termOrder, bool noMarking)
{
  Polynomial lift(p.getRing());
  Polynomial r(p.getRing());

  Monomial marked=p.getMarked().m;

  for(PolynomialSet::iterator i=l.begin();i!=l.end();i++)
    i->scaleMarkedCoefficientToOne();

  while(!p.isZero())
    {
      p.mark(termOrder);

      /*      fprintf(Stderr,"Polynomial:\n");
      AsciiPrinter(Stderr).printPolynomial(p);
      fprintf(Stderr,"\n");
      fprintf(Stderr,"Remainder:\n");
      AsciiPrinter(Stderr).printPolynomial(r);
      fprintf(Stderr,"\n");
      */

      Term initial=p.getMarked();

      PolynomialSet::const_iterator i;
      PolynomialSet::const_iterator iLift=lLift.begin();

      for(i=l.begin();i!=l.end();i++)
        {
          if(i->getMarked().m.exponent.divides(initial.m.exponent))break;
          iLift++;
        }
      if(i!=l.end())
        {
          Term s(initial.c,Monomial(p.getRing(),initial.m.exponent-i->getMarked().m.exponent));
          p-=((*i)*s);
          //          lift+=(*iLift)*s;
          lift.madd(s,*iLift);
        }
      else
        {
          p-=initial;
          r+=initial;
        }
    }
  if(!noMarking)lift.mark(marked);

  return lift;
}
#else

Polynomial divisionLift(Polynomial p, PolynomialSet l, PolynomialSet lLift, TermOrder const &termOrder, bool noMarking)
{
  Polynomial lift(p.getRing());
  Polynomial r(p.getRing());

  Monomial marked=p.getMarked().m;

  for(PolynomialSet::iterator i=l.begin();i!=l.end();i++)
    i->scaleMarkedCoefficientToOne();

  PolynomialSet coefficientPolynomials(p.getRing());
  for(PolynomialSet::iterator i=l.begin();i!=l.end();i++)coefficientPolynomials.push_back(p.getRing());


  while(!p.isZero())
    {
      p.mark(termOrder);

      /*      fprintf(Stderr,"Polynomial:\n");
      AsciiPrinter(Stderr).printPolynomial(p);
      fprintf(Stderr,"\n");
      fprintf(Stderr,"Remainder:\n");
      AsciiPrinter(Stderr).printPolynomial(r);
      fprintf(Stderr,"\n");
      */

      Term initial=p.getMarked();

      PolynomialSet::const_iterator i;
      PolynomialSet::const_iterator iLift=lLift.begin();
      PolynomialSet::iterator iCoeff=coefficientPolynomials.begin();

      for(i=l.begin();i!=l.end();i++)
        {
          if(i->getMarked().m.exponent.divides(initial.m.exponent))break;
          iLift++;
          iCoeff++;
        }
      if(i!=l.end())
        {
          Term s(initial.c,Monomial(p.getRing(),initial.m.exponent-i->getMarked().m.exponent));
          p-=((*i)*s);
          //          lift+=(*iLift)*s;
//          lift.madd(s,*iLift);
          *iCoeff+=s;
        }
      else
        {
          p-=initial;
          r+=initial;
        }
    }

  PolynomialSet::const_iterator iCoeff=coefficientPolynomials.begin();
  for(PolynomialSet::const_iterator i=lLift.begin();i!=lLift.end();i++)
    {
      lift+=*i* *iCoeff;
      iCoeff++;
    }


  if(!noMarking)lift.mark(marked);

  return lift;
}
#endif

bool isIdealContainedInIdeal(PolynomialSet const &generators, PolynomialSet const &groebnerBasis)
{
  for(PolynomialSet::const_iterator i=generators.begin();i!=generators.end();i++)
    {
      //      if(!division(*i,groebnerBasis,LexicographicTermOrder()).isZero())return false;
      if(!division(*i,groebnerBasis,StandardGradedLexicographicTermOrder()).isZero())return false;
      //      fprintf(Stderr,".\n");
    }
  return true;
}

bool areIdealsEqual(PolynomialSet const &groebnerBasis1, PolynomialSet const &groebnerBasis2)
{
  return isIdealContainedInIdeal(groebnerBasis1,groebnerBasis2)
    && isIdealContainedInIdeal(groebnerBasis2,groebnerBasis1);
}
