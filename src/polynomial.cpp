#include "polynomial.h"

#include <sstream>
#include <iostream>

#include "division.h"
#include "printer.h"
#include "timer.h"
#include "linalg.h"
#include "wallideal.h"
#include "division.h"//used for Polynomial::divides()

//static Timer polynomialTimer("Polynomial subtraction",10);
//static Timer polynomialTimer1("Polynomial monomial multiplication",1);

//-----------------------------------------
// Polynomial
//-----------------------------------------

Polynomial::Polynomial(const Term &t):
  marked(t),
  sugar(0),
  theRing(t.getRing())
{
//  terms[t.m]=t.c;
	if(!t.c.isZero())
	{
		terms.insert(std::pair<Monomial,FieldElement>(t.m,t.c));
		isMarkedBool=true;
	}
	else isMarkedBool=false;
}


Polynomial::Polynomial(PolynomialRing const &r):
  theRing(r),
  sugar(0),
  marked(r)
{
	isMarkedBool=false;
}

Polynomial Polynomial::half(bool secondHalf)const
{
  Polynomial ret(theRing);
  int numberOfTerms=terms.size();
  int splitIndex=numberOfTerms/2;
  int I=0;
  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++,I++)
    if((secondHalf)!=(I<splitIndex))
    {
      ret.terms.insert(ret.terms.end(),*i);
    }
  return ret;
}


void Polynomial::mark(TermOrder const &termOrder)
{
if(1)
	{
//	assert(!terms.empty());
	if(terms.empty()){isMarkedBool=true;return;}
	TermMap::const_iterator best=terms.begin();
	TermMap::const_iterator i=best;i++;
	while(i!=terms.end())
	{
		if(termOrder(best->first.exponent,i->first.exponent))best=i;
		i++;
	}
	marked=Term(best->second,best->first);
	isMarkedBool=true;
	}
else
{
	TermMap::iterator i=terms.begin();

  if(i!=terms.end())
    {
      Term best=Term(i->second,i->first);

      for(;i!=terms.end();i++)
	if(termOrder(best.m.exponent,i->first.exponent))best=Term(i->second,i->first);
      marked=best;
    }
  isMarkedBool=true;
}
}


void Polynomial::mark(Monomial const &monomial)
{
  assert(!terms.empty());
  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    if(i->first.exponent==monomial.exponent)
      {
	marked=Term(i->second,monomial);
	isMarkedBool=true;
	return;
      }
  fprintf(Stderr,"Monomial ");
  AsciiPrinter(Stderr).printMonomial(monomial);
  fprintf(Stderr," not found in ");
  AsciiPrinter(Stderr).printPolynomial(*this);
  fprintf(Stderr,"\n");
  assert(0);
}


void Polynomial::copyMarking(Polynomial const &p)
{
  mark(p.marked.m);
}


bool Polynomial::checkMarking(TermOrder const &termOrder)const
{
	if(!isMarked())return false;
	for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      //      if(!termOrder(i->first.exponent,marked.m.exponent))return false;
      if(termOrder(marked.m.exponent,i->first.exponent))return false;
    }
  return true;
}


bool Polynomial::isHomogeneous(IntegerVector const &v)const
{
  if(isZero())return true;
  return degree(v)==-degree(-v);
}


bool Polynomial::isHomogeneous()const
{
  if(isZero())return true;
  IntegerVector v=IntegerVector::allOnes(theRing.getNumberOfVariables());
  return degree(v)==-degree(-v);
}


void Polynomial::scaleMarkedCoefficientToOne()
{
	assert(isMarked());
  FieldElement a=marked.c.inverse();
  *this*=a;
  marked.c=marked.c.one();
}


int Polynomial::getNumberOfVariables()const
{
  TermMap::const_iterator i=terms.begin();
  if(i==terms.end())return 0;

  return i->first.exponent.size();
}


void Polynomial::changeNumberOfVariables(PolynomialRing const &newRing)
{
  Polynomial q(newRing);

  for(TermMap::iterator i=terms.begin();i!=terms.end();i++)
    {
      IntegerVector v=i->first.exponent;
      v.resize(newRing.getNumberOfVariables());
      FieldElement e=i->second;
      q+=Term(e,Monomial(newRing,v));
    }
  if(isMarked())marked.m.exponent.resize(newRing.getNumberOfVariables());
  terms=q.terms;
  this->theRing=newRing;//This is bad style
}


void Polynomial::madd(const Term &m, const Polynomial &p)
{
  if(p.terms.empty())return; //added May 7 2005
  int sugar2=p.getSugar()+m.m.exponent.sum();
  if(sugar2>sugar)sugar=sugar2;
  TermMap::iterator i=terms.lower_bound(Monomial(theRing,p.terms.begin()->first.exponent+m.m.exponent));

  for(TermMap::const_iterator j=p.terms.begin();j!=p.terms.end();j++)
    {
      while(i!=terms.end() && TermMapCompare()(i->first,Monomial(theRing,j->first.exponent+m.m.exponent)))i++;
      if(i==terms.end())
	{
	  terms.insert(i,TermMap::value_type(Monomial(theRing,j->first.exponent+m.m.exponent),j->second*m.c));
	}
      else
	{
	  if(!TermMapCompare()(Monomial(theRing,j->first.exponent+m.m.exponent),i->first))
	    { // they must be equal
	      FieldElement c=i->second+j->second*m.c;
	      if(c.isZero())
		{
		  TermMap::iterator oldI=i;
		  i++;
		  terms.erase(oldI);
		}
	      else
		{
		  i->second=c;
		}
	    }
	  else
	    {
	      terms.insert(i,TermMap::value_type(Monomial(theRing,j->first.exponent+m.m.exponent),j->second*m.c));
	    }
	}
    }
}

void Polynomial::operator+=(const Polynomial &p)
{
  if(p.terms.empty())return; //added May 7 2005
  if(p.getSugar()>sugar)sugar=p.getSugar();
  // fast addition
  //  TimerScope ts(&polynomialTimer);

  TermMap::iterator i=terms.lower_bound(p.terms.begin()->first);

  for(TermMap::const_iterator j=p.terms.begin();j!=p.terms.end();j++)
    {
      while(i!=terms.end() && (terms.value_comp()(*i,*j)))i++;
      if(i==terms.end())
	{
	  terms.insert(i,TermMap::value_type(j->first,j->second));
	}
      else
	{
	  if(!terms.value_comp()(*j,*i))
	    { // they must be equal
	      FieldElement c=i->second+j->second;
	      if(c.isZero())
		{
		  TermMap::iterator oldI=i;
		  i++;
		  terms.erase(oldI);
		}
	      else
		{
		  i->second=c;
		}
	    }
	  else
	    {
	      terms.insert(i,TermMap::value_type(j->first,j->second));
	    }
	}
    }
  // slow addition
  /*  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    {
      if(terms.count(i->first)==1)
        {
          terms[i->first]=terms[i->first]+i->second;
          if(terms[i->first].isZero())terms.erase(i->first);
        }
      else
        terms[i->first]=i->second;
    }
  */
}

void Polynomial::operator-=(const Polynomial &p)
{
  if(p.terms.empty())return; //added May 7 2005
  if(p.getSugar()>sugar)sugar=p.getSugar();
  // fast subtraction
  //TimerScope ts(&polynomialTimer);

  TermMap::iterator i=terms.lower_bound(p.terms.begin()->first);

  for(TermMap::const_iterator j=p.terms.begin();j!=p.terms.end();j++)
    {
      while(i!=terms.end() && (terms.value_comp()(*i,*j)))i++;
      if(i==terms.end())
	{
	  terms.insert(i,TermMap::value_type(j->first,-j->second));
	}
      else
	{
	  if(!terms.value_comp()(*j,*i))
	    { // they must be equal
	      FieldElement c=i->second-j->second;
	      if(c.isZero())
		{
		  //		  if(i->fist.exponent==marked.exponent)marked=Term);
		  TermMap::iterator oldI=i;
		  i++;
		  terms.erase(oldI);
		}
	      else
		{
		  i->second=c;
		}
	    }
	  else
	    {
	      terms.insert(i,TermMap::value_type(j->first,-j->second));
	    }
	}
    }

  // slow subtraction
  /*  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    {
      if(terms.count(i->first)==1)
        {
          terms[i->first]=terms[i->first]-i->second;
          if(terms[i->first].isZero())terms.erase(i->first);
        }
      else
        terms[i->first]=-i->second;
    }
  */
}


Polynomial operator+(const Polynomial &p, const Polynomial &q)
{
  Polynomial r(p);
  r+=q;
  return r;
}


Polynomial operator-(const Polynomial &p, const Polynomial &q)
{
  Polynomial r(p);
  r-=q;
  return r;
}


void Polynomial::operator*=(const Term &t)
{
  sugar+=t.m.exponent.sum();
#if 0
  Polynomial p(theRing);
  for(TermMap::iterator i=terms.begin();i!=terms.end();i++)
    {
      FieldElement prod=i->second;
      prod*=t.c;
      p.terms.insert(p.terms.end(),TermMap::value_type(Monomial(theRing,i->first.exponent+t.m.exponent),prod));
    }
  terms=p.terms;
#else //This might violate the C++ standard but will work for any sane implementation of STL.
  for(TermMap::iterator i=terms.begin();i!=terms.end();i++)
	  {
	  	  (IntegerVector &)(i->first.exponent)+=t.m.exponent;
	  	  i->second*=t.c;
	  }
#endif
}


void Polynomial::operator*=(const Monomial &m)
{
  sugar+=m.exponent.sum();
#if 0
  Polynomial p(theRing);
  for(TermMap::iterator i=terms.begin();i!=terms.end();i++)
    {
      FieldElement prod=i->second;
      p.terms.insert(p.terms.end(),TermMap::value_type(Monomial(theRing,i->first.exponent+m.exponent),prod));
    }
  terms=p.terms;
#else //This might violate the C++ standard but will work for any sane implementation of STL.
  for(TermMap::iterator i=terms.begin();i!=terms.end();i++)(IntegerVector &)(i->first.exponent)+=m.exponent;
#endif
}


void Polynomial::operator*=(FieldElement const &c)
{
  for(TermMap::iterator i=terms.begin();i!=terms.end();i++)
    {
      i->second*=c;
    }
}

void Polynomial::operator*=(Polynomial const &p)
{
  Polynomial r(theRing);

#if 0
  if(terms.size()>p.terms.size())
    {Polynomial q=p;q*=*this;*this=q;}
  else if(isZero())
    *this=Polynomial(theRing);
  else if(terms.size()==1)
    *this=p*Term(terms.begin()->second,terms.begin()->first);
  else
    {
      Polynomial A=half(true);
      Polynomial B=half(false);
      A*=p;
      B*=p;
      *this=A+B;
    }
#else
  if(terms.size()>p.terms.size())
    {
      Polynomial q=p;q*=*this;*this=q;
    }
  else
  {
  for(TermMap::iterator i=terms.begin();i!=terms.end();i++)
    {
      r+=p*Term(i->second,i->first);
    }
  *this=r;
  }
#endif
  }

Polynomial operator*(const Polynomial &p, Term const &t)
{
  Polynomial r(p);
  r*=t;

  return r;
}


Polynomial operator*(const Polynomial &p, Monomial const &m)
{
  Polynomial r(p);
  r*=m;

  return r;
}


Polynomial operator*(const Polynomial &p, FieldElement const &c)
{
  Polynomial r(p);
  r*=c;

  return r;
}


Polynomial operator*(const Polynomial &p, const Polynomial &q)
{
  Polynomial r(p);
  r*=q;

  return r;
}


bool Polynomial::divides(const Polynomial &p, Polynomial *result)const
{
//	debug<<p<<" BY "<<*this<<"\n";
//	{static int a;a++;assert(a<10000);}
	if(isZero())return false;
	StandardGradedLexicographicTermOrder T;
	PolynomialSet l(getRing());
	l.push_back(*this);
	l.back().mark(T);
	Polynomial P=p;
	P.mark(T);
	if(result)
	{
		PolynomialSet res(getRing());
		Polynomial r=division(P,l,T,&res);
		*result=res.front();
		return r.isZero();
	}
	Polynomial r=division(P,l,T);
	return r.isZero();
}


bool Polynomial::divides2(const Polynomial &p, Polynomial *result)const
{
//	debug<<p<<" B2Y "<<*this<<"\n";
	if(isZero())return false;
	StandardGradedLexicographicTermOrder T;
	PolynomialSet l(getRing());
	l.push_back(*this);
	l.back().mark(T);
	Polynomial P=p;
	P.mark(T);
	if(result)
	{
		PolynomialSet res(getRing());
		Polynomial r=division(P,l,T,&res);
		*result=res.front();
		return r.isZero();
	}
	Polynomial r=division(P,l,T);
	return r.isZero();
}


Polynomial Polynomial::exactlyDividedBy(const Polynomial &p)const
{
//	debug<<p<<" EXACT BY "<<*this<<"\n";
	assert(!p.isZero());

	StandardGradedLexicographicTermOrder T;
	PolynomialSet l(getRing());
	l.push_back(p);
	l.back().mark(T);
	Polynomial P=*this;
	P.mark(T);
	PolynomialSet res(getRing());
	Polynomial r=division(P,l,T,&res);
	assert(r.isZero());
	return res.front();
}


bool Polynomial::isZero()const
{
  return terms.begin()==terms.end();
}

bool Polynomial::isOne()const
{
  return terms.size()==1 && terms.begin()->first.exponent.isZero() && terms.begin()->second.isOne();
}

int Polynomial::numberOfTerms()const
{
  return terms.size();
}

IntegerVector Polynomial::exponentsSum()const
{
  IntegerVector sum(numberOfVariablesInRing());

  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    sum+=i->first.exponent;

  return sum;
}


IntegerVectorList Polynomial::exponents()const
{
  IntegerVectorList ret;

  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    ret.push_back(i->first.exponent);

  return ret;
}


IntegerVector Polynomial::greatestCommonMonomialDivisor()const
{
  assert(!isZero());
  IntegerVector ret=terms.begin()->first.exponent;
  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    ret=min(ret,i->first.exponent);
  return ret;
}


IntegerVector Polynomial::degreeVector()const
{
  IntegerVector ret(this->getRing().getNumberOfVariables());
  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    ret=max(ret,i->first.exponent);
  return ret;
}


bool Polynomial::isMonomial()const
{
  return terms.size()==1;//could it be faster to compare begin and end iterators?
}

int Polynomial::totalDegree()const
{
  int d=-1;
  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    if(i->first.exponent.sum()>d)d=i->first.exponent.sum();

  return d;
}


int64 Polynomial::degree(IntegerVector const &w)const
{
  bool first=true;
  int64 d=0;

  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      if(first||dotLong(i->first.exponent,w)>d)d=dotLong(i->first.exponent,w);
      first=false;
    }
  assert(!first);

  return d;
}


Polynomial Polynomial::homogenization(PolynomialRing const &newRing, IntegerVector const *w)const   //does not compute sugar
{
  int degree;
  Polynomial ret(newRing);

  if(w)
    degree=this->degree(*w);
  else
    degree=totalDegree();

  IntegerVector m;

  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      IntegerVector v=i->first.exponent;
      int d;
      if(w)
	d=dot(v,*w);
      else
	d=v.sum();
      IntegerVector a(v.size());

      a.grow(a.size()+1);
      v.grow(v.size()+1);
      a[a.size()-1]=degree-d;

      v+=a;
      ret+=Term(i->second,Monomial(newRing,v));
      if(isMarked())if(marked.m.exponent==i->first.exponent)m=v;
    }

  if(isMarked())
    {
	  IntegerVector v=m;
/*	  if(m.size()+1==newRing.getNumberOfVariables())
	  {
	  v.grow(m.size()+1);
	  v[m.size()]=degree-v.sum();
	  ret.mark(Monomial(newRing,v));
	  }*/
	  if(m.size()==newRing.getNumberOfVariables())
	  {
	  ret.mark(Monomial(newRing,v));
	  }
	  //      assert(m.size()==newRing.getNumberOfVariables());
//      ret.mark(Monomial(newRing,m));
	  assert(ret.isMarked());

    }
  return ret;
}


Polynomial Polynomial::torusAct(FieldVector const &w)const
{
  Polynomial ret=*this;
  int n=theRing.getNumberOfVariables();
  assert(w.size()==n);
  for(TermMap::iterator i=ret.terms.begin();i!=ret.terms.end();i++)
    {
      FieldElement c=theRing.getField().zHomomorphism(1);
      for(int j=0;j<n;j++)
	for(int k=0;k<i->first.exponent[j];k++)
	  c*=w[j];
      i->second=i->second*c;
    }

  if(isMarked())
    {
      ret.mark(marked.m);
    }
  return ret;
}


/*Polynomial Polynomial::derivative()const
{
  Polynomial ret(theRing);

  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      IntegerVector v=i->first.exponent;
      assert(v.size()==1);
      if(v[0]!=0)
	{
	  v[0]--;
	  ret+=Term(i->second*theRing.getField().zHomomorphism(v[0]+1),Monomial(theRing,v));
	}
    }

  return ret;
}*/


Polynomial Polynomial::derivative(int j)const
{
  Polynomial ret(theRing);

  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      IntegerVector v=i->first.exponent;
      if(v[j]!=0)
      {
    	  v[j]--;
    	  ret+=Term(i->second*theRing.getField().zHomomorphism(v[j]+1),Monomial(theRing,v));
      }
    }

  return ret;
}


Polynomial Polynomial::deHomogenization(PolynomialRing const &r)const
{
  Polynomial ret=*this;
  int n=numberOfVariablesInRing();
  assert(n>0);
  ret.changeNumberOfVariables(r);
  return ret;
}


Polynomial Polynomial::deHomogenizationInSameRing()const
{
  Polynomial ret(getRing());
  int n=getRing().getNumberOfVariables();
  assert(n>0);
  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      IntegerVector v=i->first.exponent;
      v[n-1]=0;
      FieldElement e=i->second;
      ret+=Term(e,Monomial(theRing,v));
    }
  if(isMarked())
  {
	  ret.marked.m.exponent=marked.m.exponent.subvector(0,n-1);
  ret.isMarkedBool=true;
  }
	  return ret;
}


int Polynomial::numberOfVariablesInRing()const
{
  return theRing.getNumberOfVariables();
//  assert(terms.size()!=0);
//  return terms.begin()->first.exponent.size();
}


void Polynomial::saturate(int variableNum)//does not compute sugar
{ // Should use greatestCommonMonomialDivisor() and operator* instead.
  if(!terms.empty())
    {
      IntegerVector smallest=terms.begin()->first.exponent;

      for(TermMap::iterator i=terms.begin();i!=terms.end();i++)
	smallest=min(smallest,i->first.exponent);

      if(variableNum!=-1)
	{
	  for(int j=0;j<smallest.size();j++)if(j!=variableNum)smallest[j]=0;
	}

      if(!smallest.isZero())
      {
      Polynomial p(theRing);
      for(TermMap::iterator i=terms.begin();i!=terms.end();i++)
	p.terms.insert(p.terms.end(),TermMap::value_type(Monomial(theRing,i->first.exponent-smallest),i->second));
      terms=p.terms;

    //for(TermMap::iterator i=terms.begin();i!=terms.end();i++)i->first.exponent-=smallest; Would be faster but does not compile



      if(isMarked())
	{
	  marked.m.exponent-=smallest;
	}
      }
    }
}


void Polynomial::computeInitialSugar()
{
  sugar=totalDegree();
}

int Polynomial::getSugar()const
{
  return sugar;
}


bool Polynomial::isMarked()const
{
	return isMarkedBool;
//  return marked.m.exponent.size()!=0;
}


bool Polynomial::isValid(int numberOfVariables)const
{
  if(!terms.empty())
    {
      if(numberOfVariables==-1)numberOfVariables=numberOfVariablesInRing();
      for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
	{
	  if(i->first.exponent.size()!=numberOfVariables)
	    {
	      fprintf(Stderr,"Polynomial::isValid failed!!!!\n");
	      return false;
	    }
	}
      if(isMarked())
	{
	  assert(marked.m.exponent.size()==numberOfVariables);
	}
    }
  return true;
}


FieldElement Polynomial::evaluate(const FieldElement &x)const
{
  FieldElement r=x-x;
  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      IntegerVector v=i->first.exponent;
      assert(v.size()==1);
      FieldElement s=theRing.getField().zHomomorphism(1);
      for(int j=0;j<v[0];j++)
	{
	  s=s*x;
	}
      r=r+i->second*s;
    }
  return r;
}

int Polynomial::maximalIndexOfVariableInSupport()const
{
  int ret=-1;

  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      IntegerVector v=i->first.exponent;
      for(int j=ret+1;j<v.size();j++)
	if(v[j])ret=j;
    }
  return ret;
}


IntegerVector Polynomial::usedVariables()const
{
	int n=theRing.getNumberOfVariables();
	IntegerVector ret(n);
	for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
	{
		ret=max(ret,i->first.exponent.supportAsZeroOneVector());
	}
	return ret;
}


Polynomial Polynomial::embeddedInto(PolynomialRing const &r2, list<int> const *chosenVariables)const
{
  Polynomial q(r2);

  if(chosenVariables)
    {
      assert(chosenVariables->size()==r2.getNumberOfVariables());
    }

  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      IntegerVector v=i->first.exponent;
      if(chosenVariables)
	v=v.subvector(*chosenVariables);
      else
	v.resize(r2.getNumberOfVariables());
      FieldElement e=i->second;
      q+=Term(e,Monomial(r2,v));
    }
  if(isMarked())
    {
      IntegerVector m=marked.m.exponent;
      //AsciiPrinter(Stderr)<<*this<<"\n";
      if(chosenVariables)
	m=m.subvector(*chosenVariables);
      else
	m.resize(r2.getNumberOfVariables());
      q.mark(Monomial(r2,m));
    }
  //  q.marked.m.exponent.resize(r2.getNumberOfVariables());
  return q;
}


Polynomial Polynomial::embeddedInto2(PolynomialRing const &r2, vector<int> const &positionOfVariables)const
{
  Polynomial q(r2);

  assert(positionOfVariables.size()==getRing().getNumberOfVariables());

  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      IntegerVector v=i->first.exponent;
      v=v.expanded(r2.getNumberOfVariables(),positionOfVariables);
      FieldElement e=i->second;
      q+=Term(e,Monomial(r2,v));
    }
  if(isMarked())
    {
      IntegerVector m=marked.m.exponent;
      m=m.expanded(r2.getNumberOfVariables(),positionOfVariables);
      q.mark(Monomial(r2,m));
    }
  return q;
}

Polynomial Polynomial::withRestrictedVariables(IntegerVector const &keepVariable, PolynomialRing const &r2)const
{
	assert(theRing.getNumberOfVariables()==keepVariable.size());
	assert(r2.getNumberOfVariables()==keepVariable.sum());

	Polynomial q(r2);
	for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
		q+=Term(i->second,Monomial(r2,i->first.exponent.subvectorSubsetBoolean(keepVariable)))	;
	return q;
}

Polynomial Polynomial::withExpandedVariables(IntegerVector const &wasKeptVariables, PolynomialRing const &r1)const
{
	assert(theRing.getNumberOfVariables()==wasKeptVariables.sum());
	assert(r1.getNumberOfVariables()==wasKeptVariables.size());

	Polynomial q(r1);
	for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
		q+=Term(i->second,Monomial(r1,i->first.exponent.expandedBoolean(wasKeptVariables)));
	return q;
}


int Polynomial::numberOfVariablesInUseConsecutive()const
{
  int ret=0;

  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      int a=i->first.exponent.indexOfLargestNonzeroEntry()+1;
      if(a>ret)ret=a;
    }

  return ret;
}


string Polynomial::toString(bool latex/*, bool mathMode*/)const
{
  stringstream s;
  /*
  if(latex && !mathMode)
    s << "$";
  */
  bool first=true;

  if(terms.empty())
    {
      s << "0";
      return s.str();
    }
  // If the polynomial has a marked term it is written first
  //   printString("_");
  IntegerVector e=getMarked().m.exponent;
  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    if(e==i->first.exponent)
      {
	s << i->second.toString(i->first.exponent.isZero(),!first,latex);
	if((!i->first.exponent.isZero())&&(!i->second.isOne())&&(!(-(i->second)).isOne()))s<<"*";
	s << i->first.toString(false,false,latex);
	first=false;
      }
  //    printString("_");
  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    if(e!=i->first.exponent)
      {
	s << i->second.toString(i->first.exponent.isZero(),!first,latex);
	if((!i->first.exponent.isZero())&&(!i->second.isOne())&&(!(-(i->second)).isOne()))s<<"*";
	s << i->first.toString(false,false,latex);
	/*	printFieldElement(i->second,i->first.exponent.isZero(),!first);
	if((!i->first.exponent.isZero())&&(!i->second.isOne())&&(!(-(i->second)).isOne()))printString("*");
	printMonomial(i->first,false,false);*/
	first=false;
      }
  /*  if(latex && !mathMode)
    s << "$";
  */
   return s.str();
}


bool Polynomial::checkExponentVectors()const
{
	int n=getRing().getNumberOfVariables();
	for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
		if(i->first.exponent.v.size()!=n)
		{
		//	log1
			{
			//	AsciiPrinter(Stderr)<<"Exponent vector length does not match ring\n"<<*this;
				assert(0);
			}
			return false;
		}
	return true;
}


double Polynomial::evaluateFloat(FloatVector const &x)const
{
	int n=getRing().getNumberOfVariables();
	assert(x.size()==n);
	double ret=0;

	for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
	{
		IntegerVector const v=i->first.exponent;
		double m=i->second.floatingPointApproximation();
		for(int i=0;i<n;i++)
		{
			for(int e=0;e<v[i];e++)m*=x[i];
		}
		ret+=m;
	}
	return ret;
}


ComplexNumber Polynomial::evaluateComplex(ComplexVector const &x)const
{
	int n=getRing().getNumberOfVariables();
	assert(x.size()==n);
	ComplexNumber ret(0,0);

	for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
	{
		IntegerVector const v=i->first.exponent;
		ComplexNumber m=i->second.floatingPointApproximation();
		for(int i=0;i<n;i++)
		{
			for(int e=0;e<v[i];e++)m*=x[i];
		}
		ret+=m;
	}
	return ret;
}


Polynomial Polynomial::withIthVariableSubstituted(PolynomialRing const &r2, int i, FieldElement const &v)const
//Polynomial Polynomial::withIthVariableSubstituted(PolynomialRing const &r2, int i, FieldElement const &v)const
{
	Polynomial ret(r2);

	for(TermMap::const_iterator j=terms.begin();j!=terms.end();j++)
	{
		Term temp=Term(j->second,Monomial(r2,j->first.exponent.withIthCoordinateRemoved(i)));
		for(int k=0;k<j->first.exponent[i];k++)
			temp*=Term(v,Monomial(r2));
		ret+=temp;
	}
	return ret;
}


Polynomial Polynomial::modularRepresentative(PolynomialRing const &r2)const
{
	Polynomial ret(r2);

	for(TermMap::const_iterator j=terms.begin();j!=terms.end();j++)
	{
		ret+=Term(j->second.modularRepresentative(r2.getField()),Monomial(r2,j->first.exponent));
	}
	return ret;
}


Polynomial Polynomial::integralRepresentative(PolynomialRing const &r2)const
{
	Polynomial ret(r2);

	for(TermMap::const_iterator j=terms.begin();j!=terms.end();j++)
	{
		ret+=Term(r2.getField().zHomomorphism(j->second.getIntegerRepresentation()),Monomial(r2,j->first.exponent));
	}
	return ret;
}


bool Polynomial::operator==(Polynomial const &q)const
{
	return !PolynomialCompare()(*this,q) && !PolynomialCompare()(*this,q);
}

//-----------------------------------------
// PolynomialSet
//-----------------------------------------

void PolynomialSet::saturate(int variableNum)
{
  for(iterator i=begin();i!=end();i++)
    i->saturate(variableNum);
}

PolynomialSet PolynomialSet::polynomialRingIntersection(PolynomialRing const &newRing, list<int> const *chosenVariables)const
{
  PolynomialSet ret(newRing);

  if(chosenVariables==0)
  {
	  for(PolynomialSet::const_iterator i=begin();i!=end();i++)
		  if(i->maximalIndexOfVariableInSupport()<newRing.getNumberOfVariables())
		  {
			  ret.push_back(i->embeddedInto(newRing));
		  }
  }
  else
  {
	  int n=theRing.getNumberOfVariables();
	  IntegerVector allowedVariables(n);
	  for(list<int>::const_iterator i=chosenVariables->begin();i!=chosenVariables->end();i++)allowedVariables[*i]=1;
	  for(PolynomialSet::const_iterator i=begin();i!=end();i++)
	  {
		  if(i->usedVariables().divides(allowedVariables))
		  {
			  ret.push_back(i->embeddedInto(newRing,chosenVariables));
		  }
	  }
  }
  return ret;
}

void PolynomialSet::changeNumberOfVariables(PolynomialRing const &r2)
{
/*  int numberOfVariables=0;

  if(n==-1)
    {
      for(PolynomialSet::const_iterator i=begin();i!=end();i++)
	{
	  if(i->getNumberOfVariables()>numberOfVariables)numberOfVariables=i->getNumberOfVariables();
	}
    }
  else
    {
      numberOfVariables=n;
    }
  //  fprintf(Stderr,"numberOFva%i\n",numberOfVariables);

  for(iterator i=begin();i!=end();i++)
    {
      i->changeNumberOfVariables(numberOfVariables);
    }
*/
	for(iterator i=begin();i!=end();i++)
	{
		i->changeNumberOfVariables(r2);
	}
	this->theRing=r2;//bad style
}


void PolynomialSet::mark(class TermOrder const &termOrder)
{
  for(iterator i=begin();i!=end();i++)
    i->mark(termOrder);
}

void PolynomialSet::markAndScale(TermOrder const &termOrder)
{
  mark(termOrder);
  scaleMarkedCoefficientsToOne();
}


void PolynomialSet::copyMarkings(PolynomialSet const &g)
{
  PolynomialSet::const_iterator j=g.begin();
  for(PolynomialSet::iterator i=begin();i!=end();i++,j++)
    i->copyMarking(*j);
}


void PolynomialSet::scaleMarkedCoefficientsToOne()
{
  for(iterator i=begin();i!=end();i++)
    i->scaleMarkedCoefficientToOne();
}


bool PolynomialSet::checkMarkings(TermOrder const &termOrder)const
{
  for(const_iterator i=begin();i!=end();i++)
    if(!i->checkMarking(termOrder))return false;

  return true;
}


bool PolynomialSet::containsInClosedGroebnerCone(IntegerVector const &v)const
{
  for(const_iterator i=begin();i!=end();i++)
    {
      if(i->degree(v)>dotLong(v,i->getMarked().m.exponent))return false;
    }
  return true;
}


bool PolynomialSet::isHomogeneous(IntegerVector const &v)const
{
  for(const_iterator i=begin();i!=end();i++)
    {
      if(!i->isHomogeneous(v))return false;
    }
  return true;
}


bool PolynomialSet::isHomogeneous()const
{
	IntegerVector v=IntegerVector::allOnes(theRing.getNumberOfVariables());
	return isHomogeneous(v);
}


void PolynomialSet::unionPolynomial(const Polynomial &p)
{
  const_iterator j;
  for(j=begin();j!=end();j++)
    if((p-*j).isZero())break;
  if(j==end())push_back(p);
}


void PolynomialSet::unionSet(const PolynomialSet &s)
{
  for(const_iterator i=s.begin();i!=s.end();i++)
      unionPolynomial(*i);
}



bool PolynomialCompare::operator()(const Polynomial &a, const Polynomial &b)const
{
  if(a.terms.size()<b.terms.size())return true;
  if(b.terms.size()<a.terms.size())return false;

  TermMap::const_iterator i=a.terms.begin();
  for(TermMap::const_iterator j=b.terms.begin();j!=b.terms.end();j++,i++)
    {
      if(LexicographicTermOrder()(i->first.exponent,j->first.exponent))return true;
      if(LexicographicTermOrder()(j->first.exponent,i->first.exponent))return false;
    }
  if((a-b).isZero())return false;
  if((a+b).isZero())return false; // we need some way of comparing field elements

/*  AsciiPrinter(Stderr).printPolynomial(a);
  fprintf(Stderr,"\n");
  AsciiPrinter(Stderr).printPolynomial(b);
*/

  i=a.terms.begin();
  for(TermMap::const_iterator j=b.terms.begin();j!=b.terms.end();j++,i++)
    {
      int s=(i->second - j->second).sign();//TODO: Make this work for fields which are not ordered
      if(s==1)return false;
      if(s==-1)return true;
    }
//  assert("Polynomial compare must be improved to handle this case"==0);
  return false;
}


PolynomialCompareMarkedTerms::PolynomialCompareMarkedTerms(TermOrder const &termOrder_):
  termOrder(termOrder_)
{
}


bool PolynomialCompareMarkedTerms::operator()(const Polynomial &a, const Polynomial &b)const
{
  return termOrder(a.getMarked().m.exponent,b.getMarked().m.exponent);
}


PolynomialCompareMarkedTermsReverse::PolynomialCompareMarkedTermsReverse(TermOrder const &termOrder_):
  termOrder(termOrder_)
{
}


bool PolynomialCompareMarkedTermsReverse::operator()(const Polynomial &a, const Polynomial &b)const
{
  return termOrder(b.getMarked().m.exponent,a.getMarked().m.exponent);
}

bool PolynomialCompareNumberOfTermsStable::operator()(const Polynomial &a, const Polynomial &b)const
{
  return a.numberOfTerms()<b.numberOfTerms();
}



int PolynomialSet::totalDegree()const
{
  int d=-1;
  for(const_iterator i=begin();i!=end();i++)
    if(d<i->totalDegree())d=i->totalDegree();

  return d;
}


void PolynomialSet::sort_()
{
  sort(PolynomialCompare());
}

void PolynomialSet::simplestPolynomialsFirst()
{
  sort(PolynomialCompareNumberOfTermsStable());
}

int PolynomialSet::numberOfVariablesInRing()const
{
  return theRing.getNumberOfVariables();
  //  assert(size()!=0);
  //  return begin()->numberOfVariablesInRing();
}


bool operator==(PolynomialSet const &a, PolynomialSet const &b)
{
  return b.isEqualTo(a);
}


bool PolynomialSet::isEqualTo(PolynomialSet const &a)const
{
  if(a.size()!=size())return false;

  PolynomialSet::const_iterator j=begin();

  for(PolynomialSet::const_iterator i=a.begin();i!=a.end();i++)
    {
      if(!(*i-*j).isZero())return false;
      j++;
    }
  return true;
}


bool PolynomialSet::isUnitIdeal()const
{
  return division(theRing.one(),*this,StandardGradedLexicographicTermOrder()).isZero();
}


IntegerVector PolynomialSet::exponentsSum()const
{
  IntegerVector sum(numberOfVariablesInRing());

  for(const_iterator i=begin();i!=end();i++)
    sum+=i->exponentsSum();

  return sum;
}


PolynomialSet PolynomialSet::markedTermIdeal()const
{
  PolynomialSet LT(theRing);

  for(const_iterator i=begin();i!=end();i++)
    {
      LT.push_back(Polynomial(i->getMarked()));
    }

  return LT;
}


void PolynomialSet::computeInitialSugar()
{
  for(iterator i=begin();i!=end();i++)
    i->computeInitialSugar();
}


bool PolynomialSet::isMarked()const
{
  for(const_iterator i=begin();i!=end();i++)
    if(!i->isMarked())return false;

  return true;
}


PolynomialSet PolynomialSet::torusAct(FieldVector const &w)const
{
  PolynomialSet ret(theRing);
  for(const_iterator i=begin();i!=end();i++)
    ret.push_back(i->torusAct(w));

  return ret;
}


PolynomialSet PolynomialSet::homogenization(PolynomialRing const &newRing, IntegerVector const *w)const
{
  PolynomialSet ret(newRing);
  for(const_iterator i=begin();i!=end();i++)
    ret.push_back(i->homogenization(newRing,w));

  return ret;
}


list<int> PolynomialSet::multiDeHomogenizationToKeep()const
{
	 int n=getRing().getNumberOfVariables();
	 FieldMatrix m=integerMatrixToFieldMatrix(rowsToIntegerMatrix(wallInequalities(*this),n),Q);
	 m.reduce();
	 return m.pivotColumns();
}

PolynomialSet PolynomialSet::multiDeHomogenization()const
{
  int n=getRing().getNumberOfVariables();
  list<int> toKeep=multiDeHomogenizationToKeep();
  PolynomialRing R2=PolynomialRing(getRing().getField(),toKeep.size());
  return embeddedInto(R2,&toKeep);
}


PolynomialSet PolynomialSet::deHomogenization()const
{
  PolynomialSet ret(PolynomialRing(theRing.getField(),theRing.getNumberOfVariables()-1));

  for(const_iterator i=begin();i!=end();i++)
    ret.push_back(i->deHomogenization(ret.getRing()));

  return ret;
}


PolynomialSet PolynomialSet::deHomogenizationInSameRing()const
{
  PolynomialSet ret(theRing);

  for(const_iterator i=begin();i!=end();i++)
    ret.push_back(i->deHomogenizationInSameRing());

  return ret;
}


bool PolynomialSet::isValid()const
{
  if(size()!=0)
    {
      int n=numberOfVariablesInRing();
      for(const_iterator i=begin();i!=end();i++)
	if(!i->isValid(n))return false;
    }

  return true;
}


PolynomialSet PolynomialSet::embeddedInto(PolynomialRing const &r2, list<int> const *chosenVariables)const
{
  PolynomialSet ret(r2);
  for(const_iterator i=begin();i!=end();i++)
    ret.push_back(i->embeddedInto(r2,chosenVariables));

  return ret;
}


PolynomialSet PolynomialSet::embeddedInto2(PolynomialRing const &r2, vector<int> const &positionsOfVariables)const
{
  PolynomialSet ret(r2);
  for(const_iterator i=begin();i!=end();i++)
    ret.push_back(i->embeddedInto2(r2,positionsOfVariables));

  return ret;
}


PolynomialSet PolynomialSet::withRestrictedVariables(IntegerVector const &keepVariable, PolynomialRing const &r2)const
{
	assert(theRing.getNumberOfVariables()==keepVariable.size());
	assert(r2.getNumberOfVariables()==keepVariable.sum());
	PolynomialSet ret(r2);
	for(const_iterator i=begin();i!=end();i++)
		ret.push_back(i->withRestrictedVariables(keepVariable,r2));
	return ret;
}

PolynomialSet PolynomialSet::withExpandedVariables(IntegerVector const &wasKeptVariables, PolynomialRing const &r1)const
{
	assert(theRing.getNumberOfVariables()==wasKeptVariables.sum());
	assert(r1.getNumberOfVariables()==wasKeptVariables.size());
	PolynomialSet ret(r1);
	for(const_iterator i=begin();i!=end();i++)
		ret.push_back(i->withExpandedVariables(wasKeptVariables,r1));
	return ret;
}

int PolynomialSet::numberOfVariablesInUseConsecutive()const
{
  int ret=0;

  for(const_iterator i=begin();i!=end();i++)
    {
      int a=i->numberOfVariablesInUseConsecutive();
      if(a>ret)ret=a;
    }

  return ret;
}


int PolynomialSet::totalNumberOfTerms()const
{
  int ret=0;
  for(const_iterator i=begin();i!=end();i++)
    ret+=i->numberOfTerms();
  return ret;
}

list<IntegerVectorList> PolynomialSet::exponents()const
{
  list<IntegerVectorList> ret;
  for(const_iterator i=begin();i!=end();i++)ret.push_back(i->exponents());
  return ret;
}

FloatVector PolynomialSet::evaluateFloat(FloatVector const &x)const
{
	assert(x.size()==numberOfVariablesInRing());
	FloatVector ret(size());
	int I=0;
	for(const_iterator i=begin();i!=end();i++,I++)
		ret[I]=i->evaluateFloat(x);
	return ret;
}


ComplexVector PolynomialSet::evaluateComplex(ComplexVector const &x)const
{
	assert(x.size()==numberOfVariablesInRing());
	ComplexVector ret(size());
	int I=0;
	for(const_iterator i=begin();i!=end();i++,I++)
		ret[I]=i->evaluateComplex(x);
	return ret;
}


PolynomialSet PolynomialSet::withIthVariableSubstituted(PolynomialRing const &r2, int i, FieldElement const &v)const
{
	PolynomialSet ret(r2);
	for(const_iterator j=begin();j!=end();j++)
		ret.push_back(j->withIthVariableSubstituted(r2,i,v));
	return ret;
}


void PolynomialSet::removeZeros()
{
  for(PolynomialSet::iterator i=begin();i!=end();i++)
    {
      if(i->isZero())
        {
          PolynomialSet::iterator j=i;
          j++;
          erase(i);
          j--;
          i=j;
        }
    }
}


void PolynomialSet::removeDuplicates()
{
  sort_();
  unique();
}


PolynomialSet PolynomialSet::modularRepresentative(PolynomialRing const &r2)const
{
	PolynomialSet ret(r2);
	for(PolynomialSet::const_iterator i=begin();i!=end();i++)
		ret.push_back(i->modularRepresentative(r2));
	return ret;
}
