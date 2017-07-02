#include "buchberger.h"

#include <set>
#include <algorithm>
#include <iostream>

#include "division.h"
#include "printer.h"
#include "timer.h"
#include "parser.h"
#include "log.h"

#include "gebauermoeller.h"

static Timer buchbergerTimer("Buchberger",10);


Polynomial sPolynomial(Polynomial a, Polynomial b)
{
  bool comments=false;
  if(comments)
    {
      AsciiPrinter(Stderr).printString("S(");
      AsciiPrinter(Stderr).printPolynomial(a);
      AsciiPrinter(Stderr).printString(",");
      AsciiPrinter(Stderr).printPolynomial(b);
      AsciiPrinter(Stderr).printString(")=");
    }

  //marked coefficient of a and b must be one

  IntegerVector ina=a.getMarked().m.exponent;
  IntegerVector inb=b.getMarked().m.exponent;

  IntegerVector L=max(ina,inb);

  if(comments)
    AsciiPrinter(Stderr).printVector(L);


  FieldElement const &f=a.getMarked().c;


  Polynomial A=a;A*=Monomial(a.getRing(),L-ina);
  Polynomial B=b;B*=Monomial(b.getRing(),L-inb);

  if(comments)
    {
      AsciiPrinter(Stderr).printPolynomial(A-B);
      AsciiPrinter(Stderr).printString("\n");
    }

  return A-B;
}


// Simple Buchberger

void buchbergerChain(PolynomialSet *g, TermOrder const &termOrder, bool allowSaturation);
void buchberger2(PolynomialSet *g, TermOrder const &termOrder, bool allowSaturation);
void buchberger/*Simple*/(PolynomialSet *g, TermOrder const &termOrder, bool allowSaturation)
{
	int numberOfReductions;

		return buchberger2(g, termOrder, allowSaturation);
	//return buchbergerChain(g, termOrder, allowSaturation);
	PolynomialRing theRing=g->getRing();
  //  log2 fprintf(Stderr,"ENTERING buchberger\n");
  TimerScope ts(&buchbergerTimer);
  PolynomialSet sPolynomials(theRing);

  for(PolynomialSet::const_iterator i=g->begin();i!=g->end();i++)
    if(!i->isZero())sPolynomials.push_back(*i); // It is safe and useful to ignore the 0 polynomial

  if(allowSaturation)sPolynomials.saturate();
  sPolynomials.markAndScale(termOrder);

  *g=PolynomialSet(theRing);

  while(!sPolynomials.empty())
    {
	  //pout<<int(sPolynomials.size())<<"\n";///////////////////////
      Polynomial p=*sPolynomials.begin();
      sPolynomials.pop_front();

      p=division(p,*g,termOrder);
      numberOfReductions++;
      if(!p.isZero())
        {
    //	  pout<<"NONZERO\n";
    	  if(allowSaturation)p.saturate();
          p.mark(termOrder);
          p.scaleMarkedCoefficientToOne();
	  bool isMonomial=p.isMonomial();
          for(PolynomialSet::const_iterator i=g->begin();i!=g->end();i++)
	    if((!isMonomial) || (!i->isMonomial())) // 2 % speed up!
            {
              if(!relativelyPrime(i->getMarked().m.exponent,p.getMarked().m.exponent))
                {
                  Polynomial s=sPolynomial(*i,p);
                  s.mark(termOrder); // with respect to some termorder
                  s.scaleMarkedCoefficientToOne();
                  sPolynomials.push_back(s);
                }
            }

          if(1){//Mora trick (from his book). This is not exactly what Mora suggests, but it seems to work. (Book1, Part 3, Remark 3.6.1)
        	  for(PolynomialSet::iterator i=g->begin();i!=g->end();)
        		  if(p.getMarked().m.exponent.divides(i->getMarked().m.exponent))
        		  {
        			  PolynomialSet::iterator j=i;j++;
        			  Polynomial q=sPolynomial(*i,p);
        			  if(allowSaturation)q.saturate();
        			  if(!q.isZero())
        				  {
        				  sPolynomials.push_back(q);
        			  }
        			  g->erase(i);
//        			  debug<<"erasing\n";

        			  i=j;
        		  }
        		  else
        			  i++;
          }



          g->push_back(p);
	  {
	    static int t;
	    t++;
	    //	    if((t&31)==0)fprintf(Stderr," gsize %i  spolys:%i\n",g->size(),sPolynomials.size());
	  }
        }
    //  else
    //	  pout<<"ZERO\n";
    }
  //log2  fprintf(Stderr," buchberger minimize\n");
  minimize(g);
  //log2 fprintf(Stderr," buchberger autoreduce\n");
  autoReduce(g,termOrder);
  //log2 fprintf(Stderr,"LEAVING buchberger\n\n");

  cerr<<"NumberOfReductions: "<<numberOfReductions<<std::endl;
}



// Buchberger with chain criterion

struct ChainPair
{
  PolynomialSet::const_iterator a,b;
  int A,B;
  IntegerVector lcm;
  ChainPair(PolynomialSet::const_iterator const &a_,PolynomialSet::const_iterator const &b_,int A_,int B_):
    a(a_),
    b(b_),
    A(A_),
    B(B_),
    lcm(max(a_->getMarked().m.exponent,b_->getMarked().m.exponent))
  {
	  if(B<A)
	  {
		  a=b_;
		  b=a_;
		  A=B_;
		  B=A_;
	  }
  }
  bool operator<(const ChainPair & b)const
  {
	  if(b.lcm.sum()<lcm.sum())return false;
    if(lcm.sum()<b.lcm.sum())return true;
    if(b.lcm<lcm)return true;
    if(lcm<b.lcm)return false;
    if(A<b.A)return true;
    if(A>b.A)return false;
    if(B<b.B)return true;
    if(B>b.B)return false;
    return false;
//    assert(0);
  }
};

typedef set<ChainPair> ChainPairList;

static bool canBeRemoved(ChainPairList const &P, IntegerVector const &lcm, PolynomialSet::const_iterator i, PolynomialSet::const_iterator l, int I, int L)
{
	for(ChainPairList::const_iterator t=P.begin();t!=P.end();t++)
	{
		if(t->a==i && t->b!=l && t->b->getMarked().m.exponent.divides(lcm))
		{
			if(P.count(ChainPair(l,t->b,L,t->B))==1)return true;
		}
		if(t->b==i && t->a!=l && t->a->getMarked().m.exponent.divides(lcm))
		{
			if(P.count(ChainPair(l,t->a,L,t->A))==1)return true;
		}
		if(t->a==l && t->b!=i && t->b->getMarked().m.exponent.divides(lcm))
		{
			if(P.count(ChainPair(i,t->b,I,t->B))==1)return true;
		}
		if(t->b==l && t->a!=i && t->a->getMarked().m.exponent.divides(lcm))
		{
			if(P.count(ChainPair(i,t->a,I,t->A))==1)return true;
		}
	}
	return false;

	  //  return false;
  for(ChainPairList::const_iterator t=P.begin();t!=P.end();t++)
    {
      if(t->a==i && t->b!=l && t->b->getMarked().m.exponent.divides(lcm) /*||
	 t->b==i && t->a!=l && t->a->getMarked().m.exponent.divides(lcm) ||
	 t->a==l && t->b!=i && t->b->getMarked().m.exponent.divides(lcm) ||
	 t->b==l && t->a!=i && t->a->getMarked().m.exponent.divides(lcm)*/)return true;
    }
  return false;
}

void printPairs(ChainPairList const &P)
{
//  return;
  for(ChainPairList::const_iterator t=P.begin();t!=P.end();t++)
    {
      cerr<<"("<<t->A<<","<<t->B<<")[";
      AsciiPrinter(Stderr)<<t->a->getMarked()<<","<<t->b->getMarked()<<"]<"<<t->lcm <<">";
      cerr<<endl;
    }
  cerr<<endl;
}

//void buchbergerChain(PolynomialSet *g, TermOrder const &termOrder)
void buchbergerChain(PolynomialSet *g, TermOrder const &termOrder, bool allowSaturation)
{
	int nDivisions=0;
	int nZeroDivisions=0;
	int nRemovedChain=0;
	PolynomialRing theRing=g->getRing();
  TimerScope ts(&buchbergerTimer);
//cerr<<g->size();
  {
    PolynomialSet g2(theRing);
    for(PolynomialSet::const_iterator i=g->begin();i!=g->end();i++)
      if(!i->isZero())g2.push_back(*i); // It is safe and useful to ignore the 0 polynomial

    *g=g2;
    if(allowSaturation)g->saturate();
  }
  g->markAndScale(termOrder);

  ChainPairList P;//use better data structure for this
  ChainPairList Pchecked;//use better data structure for this
  int I=0;
  for(PolynomialSet::const_iterator i=g->begin();i!=g->end();i++,I++)
    {
      int J=0;
      for(PolynomialSet::const_iterator j=g->begin();j!=i;j++,J++)
	{
    	  if(relativelyPrime(i->getMarked().m.exponent,j->getMarked().m.exponent) || (i->isMonomial() && j->isMonomial()))
    		  Pchecked.insert(ChainPair(j,i,J,I));//here
    	  else
    		  P.insert(ChainPair(j,i,J,I));//here
	}
    }

  while(!P.empty())
    {
//	  pout<<*g;
  //       cerr<<"P\n";printPairs(P);cerr<<"Pchecked\n";printPairs(Pchecked);


      PolynomialSet::const_iterator i=P.begin()->a;
      PolynomialSet::const_iterator l=P.begin()->b;
      int I=P.begin()->A;
      int L=P.begin()->B;

      if(relativelyPrime(i->getMarked().m.exponent,l->getMarked().m.exponent) || (i->isMonomial() && l->isMonomial()))
	{
	  // Pchecked.push_back(P.front());
	  // P.pop_front();
	  Pchecked.insert(*P.begin());
	  P.erase(P.begin());
	}
      else
	{
	  IntegerVector lcm=max(i->getMarked().m.exponent,l->getMarked().m.exponent);

	  if(/*canBeRemoved(P,lcm,i,l) ||*/ canBeRemoved(Pchecked,lcm,i,l,I,L))
	    {
	      //cerr<<"removin"<<endl;
	      //P.pop_front(); // This might remove elements from P with a trivial S-poly
	      P.erase(P.begin()); // This might remove elements from P with a trivial S-poly
	      nRemovedChain++;
	    }
	  else
	    {
	      Polynomial p=sPolynomial(*i,*l);

	      p.mark(termOrder);
	      p=division(p,*g,termOrder);nDivisions++;
          if(allowSaturation)p.saturate();
	      // Pchecked.push_back(P.front());
	      // P.pop_front();
	      Pchecked.insert(*P.begin());
	      P.erase(P.begin());

	      if(!p.isZero())
		{
		  p.mark(termOrder);
		  p.scaleMarkedCoefficientToOne();
		  int K=g->size();
		  g->push_back(p);
		  PolynomialSet::const_iterator k=g->end();k--;
		  int I=0;
		  for(PolynomialSet::const_iterator i=g->begin();i!=k;i++,I++)
		    {
	    	  if(relativelyPrime(i->getMarked().m.exponent,k->getMarked().m.exponent) || (i->isMonomial() && k->isMonomial()))
	    		  Pchecked.insert(ChainPair(i,k,I,K));//here
	    	  else
	    		  P.insert(ChainPair(i,k,I,K));//here
		//      P.insert(ChainPair(i,k,I,K));//here
		    }
		}
	      else
	    	  nZeroDivisions++;
	    }
	}
    }
  //  AsciiPrinter(Stderr)<<*g;
  minimize(g);
  autoReduce(g,termOrder);
  cerr<<"NDIV "<<nDivisions<<" NREMOVEDCHAIN "<<nRemovedChain<<" NZERODIVISIONS "<<nZeroDivisions<<endl;
}



/*class MyPolynomialCompare
{
	TermOrder const &termOrder;
public:
	MyPolynomialCompare(TermOrder const &termOrder_):termOrder(termOrder_)
	{

	}
  bool operator()(const Polynomial &a, const Polynomial &b)const
  {
	  if(termOrder(a.getMarked().m.exponent,b.getMarked().m.exponent)<0)return true;
	  if(termOrder(a.getMarked().m.exponent,b.getMarked().m.exponent)>0)return false;
	  return PolynomialCompare(a,b);
  }
};*/


void autoReduceNonGB(PolynomialSet *g, TermOrder const &termOrder)
{
	g->markAndScale(termOrder);
//	g->sort(PolynomialCompareMarkedTerms(termOrder));
	g->sort(PolynomialCompareMarkedTerms(TotalDegreeTieBrokenTermOrder(termOrder)));
//	debug<<"Sorted:"<<*g;

  for(PolynomialSet::iterator i=g->begin();i!=g->end();)
    {
      Polynomial temp(*i);
      PolynomialSet::iterator tempIterator=i;
      tempIterator++;
      g->erase(i);
      Monomial monomial=temp.getMarked().m;

      Polynomial remainder=division(temp,*g,termOrder);
//      debug<<remainder<<"\n";
      if(!remainder.isZero())
      {
    	  g->insert(tempIterator,remainder);
      //g->insert(tempIterator,smartDivision(temp,*g,termOrder));
    	  tempIterator--;
    	  i=tempIterator;
    	  i->mark(termOrder);
//      i->mark(monomial);
    	  i++;
      }
      else
    	  i=tempIterator;
    }
}

struct polynomialOrder
{
  TermOrder const &order;
  polynomialOrder(TermOrder const &order_):order(order_){}

  bool operator()(const Polynomial &a, const Polynomial &b)
  {
    return order(a.getMarked().m.exponent, b.getMarked().m.exponent);
  }
};

/*
 * g must be marked.
 * all initial monomials in g must be different
 */
void autoReduceSorting(PolynomialSet *g, TermOrder const &termOrder)
{
	g->sort(polynomialOrder(termOrder));
	PolynomialSet temp(g->getRing());
	for(PolynomialSet::iterator i=g->begin();i!=g->end();i++)
	{
		Monomial monomial=i->getMarked().m;
		Polynomial f=division(*i,temp,termOrder);
		f.mark(monomial);
		temp.push_back(f);
	}
	*g=temp;
}

/* Mora's example from his book:
Q[v,w,z,y,x]
{v2-xz,y2-x3,yzv-x2w}
*/
void buchberger2(PolynomialSet *g, TermOrder const &termOrder, bool allowSaturation)
{
	bool printDebug=false;
	PolynomialRing theRing=g->getRing();

//	debug<<"Buchberger2 on:"<<theRing<<*g<<"\n";


	g->markAndScale(termOrder);
	autoReduceNonGB(g,termOrder);
	g->markAndScale(termOrder);
	g->sort(PolynomialCompareMarkedTermsReverse(termOrder));
//  g->computeInitialSugar();

  if(allowSaturation)
  {
	  for(PolynomialSet::iterator i=g->begin();i!=g->end();i++)
		  i->saturate();
  }

//  debug<<"Reduced NonGB:"<<theRing<<*g<<"\n";

  //  SPairContainerPair sPolynomials(theRing,termOrder);
  	  SPairContainerOptimized sPolynomials(theRing,termOrder);
  //  SPairContainerOptimizedPacked sPolynomials(theRing,termOrder,g);

  vector<Polynomial> G;
  int indexj=0;
  for(PolynomialSet::iterator j=g->begin();j!=g->end();j++,indexj++)
  {
	  G.push_back(*j);
	  sPolynomials.updatePairs(G,indexj);
	  if(printDebug)
	  {
//		  printSPairSet(sPolynomials);
		  cout<<"-----\n";
	  }
  }

  //  printSPairSet(sPolynomials);

  int numberOfCriticalPairsConsidered=0;
  int numberOfUsefulCriticalPairs=0;

  while(!sPolynomials.isEmpty())
    {
	  pair<int,int> ij=sPolynomials.popPair();
	  Polynomial p=sPolynomial(G[ij.first],G[ij.second]);
//      if(printDebug)
//    	  cout<<"Processing"<<sPolynomials.begin()->i+1<<sPolynomials.begin()->j+1<<"\n";
      p.mark(termOrder);
      p.scaleMarkedCoefficientToOne();

      numberOfCriticalPairsConsidered++;














      p=division(p,*g,termOrder);
      if(allowSaturation)p.saturate();
      if(!p.isZero())
        {
    	  p.mark(termOrder);
          p.scaleMarkedCoefficientToOne();
          g->push_back(p);
          G.push_back(p);
          numberOfUsefulCriticalPairs++;
          log2
          {
        	  static int t;
        	  if(((++t)&=31)==0)
        		  fprintf(Stderr,"gsize:%i spolys:%i n:%i\n",(int)g->size()+1,(int)sPolynomials.size(),(int)g->getRing().getNumberOfVariables());
          }
          sPolynomials.updatePairs(G,indexj);//,truncationDegree,&grading);

          indexj++;
          if(printDebug)
          {
        	  //printSPairSet(sPolynomials);
        	  cout<<"-----\n";
          }
        }
      else
    	if(printDebug)  cout<<"zero reduction\n";
    }
  *g=PolynomialSet(theRing);
  for(PolynomialVector::const_iterator i=G.begin();i!=G.end();i++)g->push_back(*i);
//  debug<<"MINIMIZING\n";
  minimize(g);
//  debug<<"AUTOREDUCING\n";
//  debug.printPolynomialSet(*g,true);
  autoReduceSorting(g,termOrder);
//  debug.printPolynomialSet(*g,true);
//  autoReduce(g,termOrder);
//  debug<<"RETURNING\n";

  //  if(printComments)
//    fprintf(Stderr,"Number of critical pairs considered(divisions) %i Useful %i\n",numberOfCriticalPairsConsidered,numberOfUsefulCriticalPairs);
}



void minimize(PolynomialSet *g)
{//CHECK THAT THIS ROUTINE WORKS IF TWO GENERATORS HAVE THE SAME INITIAL TERM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PolynomialRing theRing=g->getRing();
  PolynomialSet ret(theRing);

  g->sort(polynomialOrder(LexicographicTermOrder()));

  for(PolynomialSet::const_iterator i=g->begin();i!=g->end();i++)
    {
      bool someDivides=false;
      for(PolynomialSet::const_iterator j=ret.begin();j!=ret.end();j++)
	//      for(PolynomialSet::const_iterator j=g->begin();j!=i;j++) //changed Feb 2009
        {
          if(j->getMarked().m.exponent.divides(i->getMarked().m.exponent))
            {
              someDivides=true;
            }
        }
      if(!someDivides)
        ret.push_back(*i);
    }

  *g=ret;
}



void autoReduce(PolynomialSet *g, TermOrder const &termOrder)
{
	/**
	 * TODO: there should be two options : supplying a termorder, or not supplying a termorder. In the latter case this routine should decide if it wants to compute one.
	 */
//  WeightTermOrder termOrder2(termorderWeight(*g));//REMOVE ME ?? JAN 2009


  for(PolynomialSet::iterator i=g->begin();i!=g->end();i++)
    {
      Polynomial temp(*i);
      PolynomialSet::iterator tempIterator=i;
      tempIterator++;
      g->erase(i);
      Monomial monomial=temp.getMarked().m;
      g->insert(tempIterator,division(temp,*g,termOrder));
      //g->insert(tempIterator,smartDivision(temp,*g,termOrder));
      tempIterator--;
      i=tempIterator;
      i->mark(monomial);
    }
}


bool isMarkedGroebnerBasis(PolynomialSet const &g)
{
  int counter=0;
  for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    {
    log2 fprintf(Stderr,"%i ",counter++);
    for(PolynomialSet::const_iterator j=i;j!=g.end();j++)
      if(!relativelyPrime(i->getMarked().m.exponent,j->getMarked().m.exponent))
	{
	  Polynomial s=sPolynomial(*i,*j);
	  if(!division(s,g,LexicographicTermOrder()).isZero())
	    {
	      log3{AsciiPrinter(Stderr)<<"Spoly("<<*i<<","<<*j<<")="<<sPolynomial(*i,*j)<<" gives remainder "<< division(s,g,LexicographicTermOrder()) <<"\n";}
	      return false;
	    }
	}
    }
  return true;
}


bool isMinimal(PolynomialSet const &g)
{
	PolynomialSet temp=g.markedTermIdeal();
	minimize(&temp);
	return temp.size()==g.size();
}


bool isReduced(PolynomialSet const &g)
{
	if(!isMinimal(g))return false;
	PolynomialSet temp=g.markedTermIdeal();
	for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
		for(TermMap::const_iterator j=i->terms.begin();j!=i->terms.end();j++)
			if(!(j->first.exponent==i->getMarked().m.exponent))
				for(PolynomialSet::const_iterator k=temp.begin();k!=temp.end();k++)
					if(k->getMarked().m.exponent.divides(j->first.exponent))return false;
	return true;
}
