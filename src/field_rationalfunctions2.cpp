#include "field_rationalfunctions2.h"

#define USEPOLYNOMIALGCD 1
/* For some reason, setting USEFACTORY to 1 breaks code. For example
./gfan _groebnerfan --stdin ~/gfan/current/examples/3x3of3x4
no longer works */

#include <memory>
#include <assert.h>
#include <gmp.h>
#include <sstream>
#include <iostream>
#include <time.h>

#include "termorder.h"
#include "division.h"
#include "buchberger.h"
#include "saturation.h"
#include "printer.h"
#include "linalg.h"

#if USEPOLYNOMIALGCD
#include "polynomialgcd.h"
#endif

#include "log.h"

int FieldElementRationalFunctions2Living;

/**
 * This field is the field of _multivariate_ rational functions.
 */
class FieldElementRationalFunction2 : public FieldElementImplementation
{
	Polynomial gcd(Polynomial a, Polynomial b)
	{
		if(a.degree(IntegerVector::standardVector(1,0))<b.degree(IntegerVector::standardVector(1,0)))
			swap(a,b);
		while(!b.isZero())
		{
			LexicographicTermOrder T;

			a.mark(T);
			b.mark(T);
			PolynomialSet B(b.getRing());B.push_back(b);

			b=division(a,B,T);
			a=B.front();
		}
		return a;
	}
	bool randomizedFactorTest(PolynomialRing const &r, Polynomial const &p, Polynomial const &q)
	//returns true if it is discovered that no factor exists. False return value does not mean anything.
	{
		return false;//DISABLED
		vector<Polynomial> l;
		l.push_back(p);
		l.push_back(q);
		int n=r.getNumberOfVariables();
		for(int i=0;i<n;i++)
		{
			int numberOfRetries=10;
			retry:
			numberOfRetries--;
			if(numberOfRetries<0) return false;//<----------------------------------------
			FieldVector val(r.getField(),n-1);
			for(int j=0;j<n-1;j++)val[j]=r.getField().zHomomorphism((rand()&127)+1);
			vector<Polynomial> l2;
			PolynomialRing r2(r.getField(),1);

		for(int I=0;I<2;I++)
			{
				Polynomial temp(r2);
				for(TermMap::const_iterator j=l[I].terms.begin();j!=l[I].terms.end();j++)
				{
					FieldElement mp=j->second;
					for(int k=0;k<n-1;k++)
						for(int l=j->first.exponent[k+(k>=i)];l>0;l--)mp*=val[k];
					IntegerVector expo(1);expo[0]=j->first.exponent[i];
					temp+=Term(mp,Monomial(r2,expo));
				}
				if(temp.isZero())goto retry;
				if(temp.degree(IntegerVector::standardVector(1,0)!=l[I].degree(IntegerVector::standardVector(n,i)))||
						temp.degree(-IntegerVector::standardVector(1,0)!=l[I].degree(-IntegerVector::standardVector(n,i))))
					goto retry;
				l2.push_back(temp);
			}
			//gcd
			if(gcd(l2[0],l2[1]).numberOfTerms()!=1)
				return false;
		}
		return true;
	}


 public:
  Polynomial p,q;
  FieldElementRationalFunction2(FieldImplementation &a):
    FieldElementImplementation(a),
    p(((FieldRationalFunctions2Implementation*)&a)->getPolynomialRing()),
    q(Term(((FieldRationalFunctions2Implementation*)&a)->getPolynomialRing().getField().zHomomorphism(1),Monomial(((FieldRationalFunctions2Implementation*)&a)->getPolynomialRing())))
  {
    FieldElementRationalFunctions2Living++;
  }
  FieldElementRationalFunction2(FieldImplementation &a,int n_):
    FieldElementImplementation(a),
    p(Term(((FieldRationalFunctions2Implementation*)&a)->getPolynomialRing().getField().zHomomorphism(n_),Monomial(((FieldRationalFunctions2Implementation*)&a)->getPolynomialRing()))),
    q(Term(((FieldRationalFunctions2Implementation*)&a)->getPolynomialRing().getField().zHomomorphism(1),Monomial(((FieldRationalFunctions2Implementation*)&a)->getPolynomialRing())))
    {
      if(n_==0)p=Polynomial(((FieldRationalFunctions2Implementation*)&a)->getPolynomialRing());
      FieldElementRationalFunctions2Living++;
    }
  void normalize(bool doExpensiveGCD=true)
  {
	  if(p.isZero()){q=p.getRing().one();return;}
	  //STEP 1: Scale leading term of denominator to 1
	  LexicographicTermOrder T;
	  q.mark(T);
	  FieldElement s=q.getMarked().c;
	  if(!s.isOne())
	  {
		  s=s.inverse();
		  q*=s;
		  p*=s;
	  }

	  //STEP 2: If denominator is monomial, then everything is easy
	  if(q.isMonomial())
	  {
		  if(q.totalDegree()==0)return;
		  if(p.isMonomial())
		  {
			  p.mark(T);
			  Monomial s(p.getRing(),-min(p.getMarked().m.exponent,q.getMarked().m.exponent));
			  p*=s;
			  q*=s;
			  return;
		  }
	  }

	  //STEP 3: We factor out monomial factors in both p and q and remember them for later
	  Monomial qm(q.getRing(),q.greatestCommonMonomialDivisor());
	  q.saturate();
	  Monomial pm(p.getRing(),p.greatestCommonMonomialDivisor());
	  p.saturate();

	  IntegerVector A=min(qm.exponent,pm.exponent);
	  qm.exponent-=A;
	  pm.exponent-=A;

	  //STEP 4: If one of our polynomials is a monomial, then we can find no more common factors
	  if(q.isMonomial()||p.isMonomial())
	  {
		  p*=pm;
		  q*=qm;
		  return;
	  }

	  // STEP 5: We check special cases where p divides q or q divides p.
	  q.mark(T);
	  p.mark(T);
	  {
		  //if q divides p then we don't need to compute gcd
		  PolynomialSet R(q.getRing());
		  PolynomialSet Q(q.getRing());
		  Q.push_back(q);
		  if(division(p,Q,T,&R).isZero())
		  {
			  //pout<<p<<"/"<<q<<"="<<*R.begin()<<"\n";

			  p=*R.begin()*pm;
			  q=Term(qm.getRing().getField().zHomomorphism(-1),qm);
			  return;
		  }
	  }
	  {
		  //if p divides q then we don't need to compute gcd
		  PolynomialSet R(p.getRing());
		  PolynomialSet P(p.getRing());
		  P.push_back(p);
		  if(division(q,P,T,&R).isZero())
		  {
			  //pout<<p<<"/"<<q<<"="<<*R.begin()<<"\n";

			  q=*R.begin()*qm;
			  p=Term(qm.getRing().getField().zHomomorphism(-1),pm);
			  return;
		  }
	  }

	  // STEP 6: This should never happen....
	  if((q-p).isZero())
	  {
		  p=Term(pm.getRing().getField().zHomomorphism(1),pm);
		  q=Term(qm.getRing().getField().zHomomorphism(1),qm);
		  return;
	  }

	  // STEP 7: Try to prove that no factor exists by random substitution (currently disabled)
	  if(randomizedFactorTest(p.getRing(),p,q))
	  {
		  //		  pout<<"FACTORTESTWORKS\n";
		  p*=pm;
		  q*=qm;
		  return;
	  }
	  //	  assert(!q.isZero());
	  /*	  if(p.isZero()){
		  q=q.getRing().one();
		  return;
	  }
	   */
//	  pout<<"p"<<p<<"q"<<q<<"\n";//<<p*q<<"\n";
	  static int FACTOR,NOFACTOR;

	  // STEP 8: Finally we have to do a gcd computation somehow
	  if(!doExpensiveGCD){p*=pm;
	  q*=qm;return;}

	  #if USEPOLYNOMIALGCD
	  {//computes gcd

//		  debug<<"gcd\n";
//debug<<"BLA\n";

//		  long start=clock();

		  Polynomial r=polynomialGCD(p,q);
//		  if(clock()-start>CLOCKS_PER_SEC)
//		  {
//			  debug<<"GCD TAKES TOO LONG:\n"<<p.getRing()<<"\n{"<<p<<",\n"<<q<<"}\n"<<r<<"\n\n";
//		  }


//		  cerr<<"RETURNNNN---------------------===============\n";
//		  debug<<p<<"\n";
//		  debug<<q<<"\n";
//		  debug<<r<<"\n";
		  PolynomialSet R(p.getRing());
		  R.push_back(r);
		  R.mark(T);
		//  PolynomialSet pp(p.getRing());
		  p.mark(T);
		  q.mark(T);
//		  debug<<p<<"\n";
//		  debug<<q<<"\n";
		  PolynomialSet Q(p.getRing());
		  PolynomialSet P(p.getRing());
		  division(p,R,T,&P);
		  division(q,R,T,&Q);
		  p=P.front();
		  q=Q.front();
//		  debug<<p<<"\n";
//		  debug<<q<<"\n";
		  p*=pm;
		  q*=qm;

//		  debug<<"dcg\n";
//		  debug<<p<<"\n";
//		  debug<<q<<"\n";
//		  cerr<<"RETURNNNN\n";
	  }
	  #else
	  {//computes lcm
	  PolynomialSet a(p.getRing());
	  PolynomialSet b(q.getRing());
	  a.push_back(p);
	  b.push_back(q);
//	  cerr<<"Ideal Intersection\n";
	  PolynomialSet c=idealIntersection(a,b);
//	  debug<<c;
//	  cerr<<"Done Ideal Intersection\n";
	  if(c.size()!=1)
	  {
		  //		  pout<<a<<b<<c;
		  assert(c.size()==1);
	  }
	  //Polynomial
	  Polynomial r=*(c.begin());


	  //	  LexicographicTermOrder T;
	  {
		  PolynomialSet Q(p.getRing());
		  PolynomialSet P(p.getRing());
		  r.mark(T);
		  a.mark(T);
		  b.mark(T);
		  division(r,b,T,&P);
		  division(r,a,T,&Q);
		  if(p.numberOfTerms()!=P.begin()->numberOfTerms()){FACTOR++;//pout<<p<<"/"<<q<<":"<<*P.begin()<<"/"<<*Q.begin()<<"\n";
		  }
		  else NOFACTOR++;
		  p=*P.begin();
		  q=*Q.begin();
		  //pout<<"r="<<r<<"\n";
		  //pout<<"p="<<p<<"\n";
		  //pout<<"q="<<q<<"\n";
		  p*=pm;
		  q*=qm;
	  }
	  }
#endif
	 // if((NOFACTOR&(256*2-1))==0)pout<<"FAC"<<FACTOR<<","<<NOFACTOR<<"\n";
  }
  FieldElementRationalFunction2(FieldImplementation &a, Polynomial const &p_, Polynomial const &q_):
    FieldElementImplementation(a),
    p(p_),
    q(q_)
    {

      FieldElementRationalFunctions2Living++;
    }
  virtual ~FieldElementRationalFunction2()
    {
      FieldElementRationalFunctions2Living--;
    }
  FieldElementRationalFunction2& operator=(const FieldElementRationalFunction2& a)
    {
      assert(0);
      return *this;
    }

  void operator*=(const FieldElementImplementation &a)
    {
      const FieldElementRationalFunction2 *A=(const FieldElementRationalFunction2*)&a;
      assert(A);
#if 1
/*      debug<<"----\n";
      debug<<p<<" - "<<A->q<<"\n";
      debug<<q<<" - "<<A->p<<"\n";
      debug<<(q*A->p)<<" - "<<p*A->q<<"\n";
*/
      Polynomial fac1=polynomialGCD(p,A->q);
      Polynomial fac2=polynomialGCD(q,A->p);

//      Polynomial pTemp=p;
//      Polynomial qTemp=q;

//    	  p*=A->p;
//    	  q*=A->q;
//    	  normalize();
//#if 1
//    	  int totalTerms=p.numberOfTerms()+q.numberOfTerms();

 //   	  p=pTemp;
 //   	  q=qTemp;
#if 0
      p*=A->p;
      q*=A->q;

//      debug<<"A\n";
      bool err=false;
      err|=!fac1.divides(p,&p);
      err|=!fac1.divides(q,&q);
      err|=!fac2.divides(p,&p);
      err|=!fac2.divides(q,&q);
//      debug<<"B\n";
      assert(!err);
#else
      Polynomial tp=p.exactlyDividedBy(fac1);
      Polynomial tAq=A->q.exactlyDividedBy(fac1);
      Polynomial tq=q.exactlyDividedBy(fac2);
      Polynomial tAp=A->p.exactlyDividedBy(fac2);

      p=tp*tAp;
      q=tq*tAq;
#endif
      normalize(false);

/*      if(totalTerms!=p.numberOfTerms()+q.numberOfTerms())
      {
//    	  debug<<pTemp<<"/"<<qTemp<<"\n";
//    	  debug<<A->p<<"/"<<A->q<<"\n";
    	  p=pTemp;
    	  q=qTemp;
    	  p*=A->p;
    	  q*=A->q;
    	  normalize();
//    	  debug<<A->p<<"/"<<A->q<<"\n";

  //  	  assert(0);
      }*/
//#endif

#else
      p*=A->p;
      q*=A->q;
      normalize();
#endif
    }
  void operator+=(const FieldElementImplementation &a)
    {
      const FieldElementRationalFunction2 *A=(const FieldElementRationalFunction2*)&a;
      assert(A);

      p=p*A->q+A->p*q;
      q=A->q*q;
      normalize();
    }
  void madd(const FieldElementImplementation &a,const FieldElementImplementation &b)
    {
      const FieldElementRationalFunction2 *A=(const FieldElementRationalFunction2*)&a;
      const FieldElementRationalFunction2 *B=(const FieldElementRationalFunction2*)&b;
      assert(A);
      assert(B);
#if 0
      p=p*(A->q*B->q)+(A->p*B->p)*q;
      q=A->q*B->q*q;
      normalize();
#else
      Polynomial fac1=polynomialGCD(B->p,A->q);
      Polynomial fac2=polynomialGCD(B->q,A->p);

      Polynomial tBp=B->p.exactlyDividedBy(fac1);
      Polynomial tAq=A->q.exactlyDividedBy(fac1);
      Polynomial tBq=B->q.exactlyDividedBy(fac2);
      Polynomial tAp=A->p.exactlyDividedBy(fac2);

      p=p*(tBq*tAq)+(tBp*tAp)*q;
      q=tBq*tAq*q;
      normalize();
#endif
    }
  FieldElementRationalFunction2 *one() const;
  bool isZero() const
    {
      return p.isZero();
    }

  FieldElementRationalFunction2 *sum(const FieldElementImplementation &b)const
    {
      const FieldElementRationalFunction2 *B=(const FieldElementRationalFunction2*)&b;

      Polynomial fac1=polynomialGCD(q,B->q);
      Polynomial nq=B->q.exactlyDividedBy(fac1)*q;
      Polynomial np=p*B->q.exactlyDividedBy(fac1)+B->p*q.exactlyDividedBy(fac1);
      // It is still possible that more common factors exist. These must appear as factors of fac1
      FieldElementRationalFunction2 *r= new FieldElementRationalFunction2(*getField(),np,nq);
      r->normalize();
//      FieldElementRationalFunction2 *r= new FieldElementRationalFunction2(*getField(),p*B->q-B->p*q,B->q*q);
      return r;
#if 0
      const FieldElementRationalFunction2 *B=(const FieldElementRationalFunction2*)&b;

//      pout<<B->q<<"------"<<q<<"\n";

      Polynomial quotient(q.getRing());
      if(B->q.divides(q,&quotient))
      {
//          pout<<B->q<<"+-----"<<q<<"--"<<quotient<<"\n";
          return new FieldElementRationalFunction2(*getField(),p-B->p*quotient,q);
      }
      else if(q.divides(B->q,&quotient))
      {
          return new FieldElementRationalFunction2(*getField(),B->p-p*quotient,B->q);
      }


      FieldElementRationalFunction2 *r= new FieldElementRationalFunction2(*getField(),p*B->q+B->p*q,B->q*q);

      return r;
#endif
    }
  FieldElementRationalFunction2 *difference(const FieldElementImplementation &b)const
    {
      const FieldElementRationalFunction2 *B=(const FieldElementRationalFunction2*)&b;

      Polynomial fac1=polynomialGCD(q,B->q);
      Polynomial nq=B->q.exactlyDividedBy(fac1)*q;
      Polynomial np=p*B->q.exactlyDividedBy(fac1)-B->p*q.exactlyDividedBy(fac1);
      // It is still possible that more common factors exist. These must appear as factors of fac1
      FieldElementRationalFunction2 *r= new FieldElementRationalFunction2(*getField(),np,nq);
      r->normalize();
//      FieldElementRationalFunction2 *r= new FieldElementRationalFunction2(*getField(),p*B->q-B->p*q,B->q*q);
      return r;
    }
  FieldElementRationalFunction2 *negation()const
    {
      FieldElementRationalFunction2 *r= new FieldElementRationalFunction2(*getField(),p-p-p,q);

      return r;
    }
  FieldElementImplementation *inverse()const
  {
    if(isZero())
      {
	AsciiPrinter P(Stderr);
	P.printString("Error inverting FieldElement: ");
	P.printPolynomial(p);
	P.printString(" ");
	P.printPolynomial(q);
	//	P.printFieldElement(*this);
	P.printString("\n");
	assert(0);
      }

    FieldElementRationalFunction2 *r= new FieldElementRationalFunction2(*getField(),q,p);

    return r;
  }

  int sign()const
  {
	  assert(0);//not an ordered field (yet)
    if(isZero())return 0;
    return p.terms.rbegin()->second.sign();
  }

  static string LaTeXTranslator(const string &s)
  {
	  assert(0);//not supported yet
/*    int startIndex=0;
    string sign;
    if(s[0]=='-')
      {
	sign=string("-");
	startIndex=1;
      }
    int slashIndex=-1;
    for(int i=startIndex;i<s.length();i++)if(s[i]=='/')slashIndex=i;
    if(slashIndex==-1)
      return string(s);

    return sign+string("{").append(s,startIndex,slashIndex-startIndex)+string("\\over ").append(s,slashIndex+1,s.length()-slashIndex-1)+string("}");
*/
	  }

  std::string toString(bool writeIfOne=true, bool alwaysWriteSign=false, bool latexMode=false /*, bool mathMode=true*/) const
  {
    stringstream s;

    //    s << "(" <<p.toString(latexMode) << "/" << q.toString(latexMode) << ")";
    s << "(";
    s<<p.toString(latexMode);
    s << "/";
    s << q.toString(latexMode);
    s << ")";
    return s.str();
  }

  FieldElementRationalFunction2 *copy()const
  {
    FieldElementRationalFunction2 *r= new FieldElementRationalFunction2(*getField());
    r->p=p;
    r->q=q;

    return r;
  }
};

PolynomialRing FieldRationalFunctions2Implementation::getPolynomialRing()const
{
  return thePolynomialRing;
}

bool FieldRationalFunctions2Implementation::isRationals()const
{
  return false;
}

int FieldRationalFunctions2Implementation::getCharacteristic()const
{
	return this->getPolynomialRing().getField().getCharacteristic();
}


FieldRationalFunctions2Implementation::FieldRationalFunctions2Implementation(PolynomialRing const &r):
  thePolynomialRing(r)
{
}

std::string FieldRationalFunctions2Implementation::toString()const
{
  stringstream s;
  s<< thePolynomialRing.getField().toString() << "("<<thePolynomialRing.toStringVariableNames()<< ")";
  return s.str();
}

FieldElementImplementation *FieldRationalFunctions2Implementation::zHomomorphismImplementation(int n)
{
  FieldElementImplementation *ret=new FieldElementRationalFunction2(*this,n);
  return ret;
}

FieldElement FieldRationalFunctions2Implementation::zHomomorphism(int n)
{
  return FieldElement(zHomomorphismImplementation(n));
}

const char *FieldRationalFunctions2Implementation::name()
{
  return "Rational functions in n variables";
}

FieldElementRationalFunction2 *FieldElementRationalFunction2::one() const
{
  return new FieldElementRationalFunction2(*getField(),1);
}



FieldRationalFunctions2::FieldRationalFunctions2(PolynomialRing const &r):
  Field(new FieldRationalFunctions2Implementation(r))
{
}


FieldElement FieldRationalFunctions2::polynomialToFraction(Polynomial const &p)
{
#if USESHAREDPTR
	FieldRationalFunctions2Implementation *imp=dynamic_cast<FieldRationalFunctions2Implementation*>(implementingObject.get());
#else
	FieldRationalFunctions2Implementation *imp=dynamic_cast<FieldRationalFunctions2Implementation*>(implementingObject);
#endif
  Polynomial q=Term(imp->getPolynomialRing().getField().zHomomorphism(1),Monomial(imp->getPolynomialRing()));

  return new FieldElementRationalFunction2(*imp, p, q);
}

/*****************************************************
 * Conversion functions
 *****************************************************/
PolynomialRing makeVariablesParameters(PolynomialRing const &r, int numberOfParameters)
{
	assert(numberOfParameters>=0);
	assert(numberOfParameters<=r.getNumberOfVariables());
	vector<string> names(numberOfParameters);
	for(int i=0;i<numberOfParameters;i++)names[i]=r.getVariableName(i);
	PolynomialRing temp(r.getField(),names);
	vector<string> names2(r.getNumberOfVariables()-numberOfParameters);
	for(int i=0;i<names2.size();i++)names2[i]=r.getVariableName(i+numberOfParameters);

	return PolynomialRing(FieldRationalFunctions2(temp),names2);
}

Polynomial makeVariablesParameters(PolynomialRing const &genericRing, Polynomial const &p)
{
	Polynomial ret(genericRing);
#if USESHAREDPTR
	FieldRationalFunctions2Implementation const *coefficientField=dynamic_cast<FieldRationalFunctions2Implementation const*>(genericRing.getField().implementingObject.get());
#else
	FieldRationalFunctions2Implementation const *coefficientField=dynamic_cast<FieldRationalFunctions2Implementation const*>(genericRing.getField().implementingObject);
#endif
	FieldRationalFunctions2 &DANGER=(FieldRationalFunctions2&)genericRing.getField();
	PolynomialRing coefRing=coefficientField->getPolynomialRing();

	for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
	{
		FieldElement c=i->second;
		IntegerVector v=i->first.exponent;
		IntegerVector coefficientExponent=v.subvector(0,p.getRing().getNumberOfVariables()-genericRing.getNumberOfVariables());
		IntegerVector monomialExponent=v.subvector(p.getRing().getNumberOfVariables()-genericRing.getNumberOfVariables(),v.size());
		FieldElement c2=DANGER.polynomialToFraction(Term( c,Monomial(coefRing, coefficientExponent)));//does the numerator not belong to a field?
		ret+=Polynomial(Term(c2,Monomial(genericRing,monomialExponent)));
	}
	return ret;
}

PolynomialSet makeVariablesParameters(PolynomialRing const &genericRing, PolynomialSet const &p)
{
	PolynomialSet ret(genericRing);
	for(PolynomialSet::const_iterator i=p.begin();i!=p.end();i++)
		ret.push_back(makeVariablesParameters(genericRing,*i));

	return ret;
}
