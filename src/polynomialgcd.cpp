/*
 * polynomialgcd.cpp
 *
 *  Created on: Apr 26, 2014
 *      Author: anders
 */


#define USEFACTORY 0

#include <iostream>
//#include <fstream>
#include <algorithm>
#include "field_rationals.h"
#include "field_zmodpz.h"

//#include "factory.h" // to be moved into a subclass
//#include "/home/anders/math/software/Singular-svn/trunk/x86_64-Linux/include/factory.h" // to be moved into a subclass

// AT THE MOMENT WE CANNOT USE FACTORY TOGETHER WITH MATRICES IN THIS FILE!!!
//#define USEFACTORY 1
#if USEFACTORY
#include "factory/factory.h"

#else
#include "linalg.h"
#include "polynomial.h"
#endif

#include "polynomialgcd.h"
#include "saturation.h"
#include "printer.h"
#include "wallideal.h"

inline int64 min(int64 a, int64 b){return a<b?a:b;}

//template std::basic_ostream<char, std::char_traits<char> >& operator<< <int>(std::basic_ostream<char, std::char_traits<char> >&, List<int> const&);

#if USEFACTORY
CanonicalForm convertToFactory(Polynomial const &a)
{
	CanonicalForm A=0;

	int n=a.getRing().getNumberOfVariables();
	for(TermMap::const_iterator i=a.terms.begin();i!=a.terms.end();i++)
	{
		CanonicalForm term;
		if(a.getRing().getField().getCharacteristic()==0)
		{
			mpq_t *c=const_cast<mpq_t*>(i->second.getGmpRationalTemporaryPointer());
			mpz_t num;
			mpz_t den;
			mpz_init_set(num,mpq_numref(*c));
			mpz_init_set(den,mpq_denref(*c));
			term=make_cf(num,den,true);
		}
		else
		{
//			term=make_cf_from_gf(i->second.getIntegerRepresentation());
			term=i->second.getIntegerRepresentation();
		}
		for(int j=0;j<n;j++)
			term*=power(Variable(j+1),i->first.exponent[j]);//HERE
		A+=term;
	}

	return A;
}

static void recursiveConvert(const CanonicalForm &f, IntegerVector &exponentVector, PolynomialRing const &r, Polynomial &result)
{
  if (f.isZero())
  {
	  return;
  }
  if ( ! f.inCoeffDomain() )
  {
    int l = f.level()-1;//HERE
    for ( CFIterator i = f; i.hasTerms(); i++ )
    {
    	exponentVector[l]=i.exp();
      recursiveConvert(i.coeff(),exponentVector,r,result);
    }
    exponentVector[l]=0;
  }
  else
  {
	  FieldElement c(r.getField());
    if ( f.isImm() )
    {
      c=r.getField().zHomomorphism(f.intval());
    }
    else
    {
    	mpz_t den;
    	gmp_denominator(f,den);
    	mpz_t num;
    	gmp_numerator(f,num);
    	FieldElement K=fieldElementFromGmpZ(&num);//Should these take the field as argument?
    	FieldElement L=fieldElementFromGmpZ(&den);
    	c=K*L.inverse();
    }
	result+=Term(c,Monomial(r,exponentVector));
  }
}


Polynomial convertFromFactory(CanonicalForm &A, PolynomialRing const &r)
{
    IntegerVector v(r.getNumberOfVariables());
    Polynomial res(r);
    recursiveConvert(A,v,r,res);
//    debug<<"RESULT:"<<res<<"\n";
    return res;
}

/*
 * The coefficients must be either rationals or Z/pZ.
 */
Polynomial polynomialGCD(Polynomial const &a, Polynomial const &b)
{
	assert(a.getRing().getField().isRationals()||a.getRing().getField().getCharacteristic()>0 /*TODO:check if field is Z/pZ */);
//	debug<<"polynomialGCD called on:\n"<<a<<"\n"<<b<<"\n";

	On(SW_USE_NTL);
    setCharacteristic( a.getRing().getField().getCharacteristic() );
    if(a.getRing().getField().isRationals())
    {
    	On( SW_USE_EZGCD );
    	On( SW_RATIONAL );
    }
    else
    {
    	On( SW_USE_EZGCD );
    	On( SW_USE_EZGCD_P );
    	Off( SW_RATIONAL );
    }

    CanonicalForm A=convertToFactory(a);
    CanonicalForm B=convertToFactory(b);
    CanonicalForm C=gcd(A,B);
    return convertFromFactory(C,a.getRing());
}
#else



/**
 * Here goes an attempt to implement multivariate gcd.
 * First we deal with the case of Z/pZ.
 * Later follows the Z case.
 */
Polynomial univariateRemainder(Polynomial a, Polynomial b)
{
//	debug<<"-------\n"<<a<<" "<<b<<"\n";
	// We assume the term ordering is global
	b.mark(LexicographicTermOrder());
	b.scaleMarkedCoefficientToOne();
	a.mark(LexicographicTermOrder());
	b.mark(LexicographicTermOrder());
	while((!a.isZero())&&(!b.isZero())&& b.getMarked().m.exponent.divides(a.getMarked().m.exponent))
	{
        Term s(a.getMarked().c,Monomial(a.getRing(),a.getMarked().m.exponent-b.getMarked().m.exponent));
//        debug<<"...."<<a<<"\n";
        a-=b*s;
		a.mark(LexicographicTermOrder());
	}


//	debug<<a<<"\n-------\n";

	return a;
}


Polynomial univariateGCD(Polynomial a, Polynomial b)
{
//	debug<<"UGCD"<<a<<" "<<b<<"=";
	if(a.totalDegree()<b.totalDegree())
	{
		swap(a,b);
	}
	while(!b.isZero())
	{
//		debug<<"\n  "<<a<<" "<<b<<"\n";
		a=univariateRemainder(a,b);

		swap(a,b);
	}
//	debug<<a<<"\n\n";
	return a;
}

PolynomialSet splitPolys2(PolynomialSet l, IntegerMatrix const &grading)
{
	PolynomialSet ret(l.getRing());

	for(PolynomialSet::const_iterator p=l.begin();p!=l.end();p++)
	{
		map<IntegerVector,TermMap> parts;

		for(TermMap::const_iterator i=p->terms.begin();i!=p->terms.end();i++)
		{
//			parts[grading.vectormultiply(i->first.exponent)][i->first]=i->second; //This line fails when not using shared pointers for implementing objects of FieldElements
			parts[grading.vectormultiply(i->first.exponent)].insert(*i);
		}

		for(map<IntegerVector,TermMap>::const_iterator i=parts.begin();i!=parts.end();i++)
		{
			Polynomial a(l.getRing());
			for(TermMap::const_iterator j=i->second.begin();j!=i->second.end();j++)
				a+=Term(j->second,Monomial(l.getRing(),j->first.exponent));
			ret.push_back(a);
		}
	}
	return ret;
}

PolynomialSet splitPoly(Polynomial const &p, IntegerVector const &gradings)//gradings is a zero-one vector
{
	PolynomialSet ret(p.getRing());
	map<IntegerVector,TermMap> parts;
	for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
	{
		parts[coordinatewiseProduct(i->first.exponent,gradings)][i->first]=i->second;//This line seems to cause a copy of a FieldElement that has only been initialised with its default constructor.
	}
	for(map<IntegerVector,TermMap>::const_iterator i=parts.begin();i!=parts.end();i++)
	{
		Polynomial a(p.getRing());
		for(TermMap::const_iterator j=i->second.begin();j!=i->second.end();j++)
			a+=Term(j->second,Monomial(p.getRing(),j->first.exponent));
		ret.push_back(a);
	}
	return ret;
}


PolynomialSet splitPolys(PolynomialSet const &s, IntegerVector const &gradings)//gradings is a zero-one vector
{
	PolynomialSet ret=s.getRing();
	for(PolynomialSet::const_iterator i=s.begin();i!=s.end();i++)
	{
		PolynomialSet temp=splitPoly(*i,gradings);
		ret.insert(ret.end(),temp.begin(),temp.end());
	}
	return ret;
}

Polynomial findFactorForRandomSubstitutionAllButI(PolynomialSet p, int i);
Polynomial findFactorForRandomSubstitution(PolynomialSet const &p, int i, FieldElement v);

static bool comparer(const Polynomial & a, const Polynomial & b)
{return a.totalDegree()> b.totalDegree();}


Polynomial withLexLargestCoefficientScaledTo(Polynomial p, FieldElement v)
{
	LexicographicTermOrder T;
	p.mark(T);
	FieldElement s=p.getMarked().c.inverse()*v;
	return p*s;
}


Polynomial naiveGCD(PolynomialSet const &p)
{
	assert(p.size());
	Polynomial ret=p.front();
	PolynomialSet::const_iterator i=p.begin();i++;
	for(;i!=p.end();i++)
	{
		PolynomialSet a(p.getRing());
		PolynomialSet b(p.getRing());
		a.push_back(ret);
		b.push_back(*i);
		PolynomialSet c=idealIntersection(a,b);
		assert(c.size()==1);
		c.front().divides(ret*(*i),&ret);
	}
	return ret;
}

bool simplifyPolynomialsForGCD(PolynomialSet const &p, PolynomialSet &newPolys, IntegerVector &support)
{
	IntegerVector dv=p.front().degreeVector();
	for(auto i=p.begin();i!=p.end();i++)dv=min(dv,i->degreeVector());
	support=dv.supportVector();

	if(support.sum()!=support.size())
	{
		int newnvars=support.sum();
		PolynomialRing r=p.getRing();
		PolynomialRing r2=r.subRing(support);
		newPolys=splitPolys(p,IntegerVector::allOnes(support.size())-support).withRestrictedVariables(support,r2);
		newPolys.saturate();
		return true;
	}
	return false;
}

PolynomialSet simplifyPolysViaHomogeneitySpace(PolynomialSet p)
{
	FieldMatrix A=spanOfHomogeneitySpaces(p);
//	debug<<A<<"\n";
	IntegerMatrix l=fieldMatrixToIntegerMatrixPrimitive(A);
//	debug<<">>"<<l<<"\n";
	return splitPolys2(p,l);
	//debug<<A<<"\n";

}




// return value 1 means that no factor was (which hopefully is equivalent to that it is unlikely that one exists).
Polynomial NonMonomialPolynomialGCDForZModP(PolynomialSet p)
{//ASSUMES THAT MONOMIAL FACTORS HAVE BEEN REMOVED
//	if(p.getRing().getNumberOfVariables()>=4)
//		debug<<"LEVEL3:"<<p.getRing()<<p<<"\n\n";
	PolynomialRing r=p.getRing();
	int n=r.getNumberOfVariables();
	assert(p.size());
	switch(r.getNumberOfVariables())
	{
	case 0:
		// If there are no variables, we return 1 (WHY?)
		return r.one();
	case 1:
	{
		PolynomialSet temp=p;
		temp.sort(comparer);
//		sort(temp.begin(),temp.end(),comparer);
/*		std::sort(
				temp.begin(), temp.end(),
				[](const Polynomial & a, const Polynomial & b) -> bool
				{return a.totalDegree()> b.totalDegree();}
		);*/
//		debug<<p<<"\n";
		auto i=p.begin();
		Polynomial ret=*i++;
		while(i!=p.end())
		{
			ret=univariateGCD(ret,*i);
			i++;
		}
		return ret;
	}
	break;
	default:
	{
#if 0
		IntegerVector dv=p.front().degreeVector();
		for(auto i=p.begin();i!=p.end();i++)dv=min(dv,i->degreeVector());
		IntegerVector support=dv.supportVector();

		if(support.sum()!=support.size())
		{// the positions where support is zero are not used in at least one polynomial.
			// Therefore it cannot we involved in the gcd.
			// The gcd must then divide each homogeneous part in all gradings w.r.t. e_i where support_i==1
//debug<<"A\n"<<support<<"\n";
			int newnvars=support.sum();
			PolynomialRing r=p.getRing();
			PolynomialRing r2=r.subRing(support);
//debug<<r2<<"\n";
//			PolynomialSet p2=splitPolys(p,support).withRestrictedVariables(support,r2);
			PolynomialSet p2=splitPolys(p,IntegerVector::allOnes(support.size())-support).withRestrictedVariables(support,r2);
//			debug<<p2<<"\n";
			p2.saturate();
//			if(p2.size()==0){debug<<p<<support<<"\n";}
			assert(p2.size());
			return NonMonomialPolynomialGCDForZModP(p2).withExpandedVariables(support,r);
		}
		// here we must pick a random substitution x_i=v
#endif
		PolynomialSet newP(p.getRing());
		IntegerVector support;
		if(simplifyPolynomialsForGCD(p,newP,support))
			return NonMonomialPolynomialGCDForZModP(newP).withExpandedVariables(support,p.getRing());
//		return naiveGCD(p);

		//		debug<<"B\n";
		for(int tries=0;tries<4;tries++)
		{
//			if(r.getNumberOfVariables()==3)debug<<"C\n";
			assert(r.getNumberOfVariables());
			// pick i
			int i=((unsigned int)rand())%(unsigned int)r.getNumberOfVariables();
			FieldElement c=r.getField().random();

//			debug<<"D\n";
			Polynomial f=findFactorForRandomSubstitution(p,i,c);
//			debug<<"E\n";

			if(f.isOne())
			{

//				debug<<p<<"\n";
//				debug<<i<<"\n";
				IntegerVector support3=IntegerVector::allOnes(r.getNumberOfVariables())-IntegerVector::standardVector(r.getNumberOfVariables(),i);
//				debug<<support3<<"HERE\n";
				//				PolynomialRing r3=r.subRing(support3);
				PolynomialSet p2=splitPolys(p,support3);//.withRestrictedVariables(support3,r3);

//				debug<<"SPLIT:"<<p2;

//				debug<<NonMonomialPolynomialGCDForZModP(p2);
				Polynomial projectionFactor=NonMonomialPolynomialGCDForZModP(p2);
				if(projectionFactor.numberOfTerms()>1)
				{
					PolynomialSet p3(r);
					for(PolynomialSet::const_iterator i=p.begin();i!=p.end();i++)
					{
//						Polynomial q(r);
//						projectionFactor.divides(*i,&q);//exact division?
						Polynomial q=i->exactlyDividedBy(projectionFactor);
						p3.push_back(q);
					}
					return projectionFactor*NonMonomialPolynomialGCDForZModP(p3);
				}

//					assert(0);
				// check if this is because the projection direction gave a factor
				// if so divide a try again. Multiply
				// if not continure , i.e. try again
				continue;
			}
			else
			{
				// we think we found a factor
//				debug<<"OTHER DIRECTION"<<findFactorForRandomSubstitutionAllButI(p,i)<<"\n";

				assert(p.size());
/*				int M=p.front().degree(IntegerVector::standardVector(n,i));
				for(PolynomialSet::const_iterator j=p.begin();j!=p.end();j++)
					M=min(M,j->degree(IntegerVector::standardVector(n,i)));
*/
				int M=findFactorForRandomSubstitutionAllButI(p,i).totalDegree();
//				M--;//////////////////// FORCE FACTOR TO CONTAIN THE ITH VARIABLE??????
//				debug<<"M:"<<M<<"\n";
				int numberOfSamples=M+2;
				vector<pair<FieldElement,Polynomial> > samples;

				int randomtries=0;
				for(int k=0;k<numberOfSamples;k++){
					FieldElement c2=r.getField().random();
					Polynomial f2=findFactorForRandomSubstitution(p,i,c2);
//					debug<<c2<<" : "<<f2<<"\n";
					if(!f2.isZero())
						samples.push_back(pair<FieldElement,Polynomial>(c2,f2));
					else
					{
						k--;
						randomtries++;
						if(randomtries>100)debug<<"TROUBLE IN GCD COMPUTATION\n";
					}
				}

//				debug<<"AA\n";
				set<IntegerVector> observedExponentsSet;
				for(vector<pair<FieldElement,Polynomial> >::const_iterator k=samples.begin();k!=samples.end();k++)
				{
					IntegerVectorList L=k->second.exponents();
//					debug.printVectorList(L);
					for(IntegerVectorList::const_iterator l=L.begin();l!=L.end();l++)
						//observedExponentsSet.insert(l->withIthCoordinateRemoved(i));
						observedExponentsSet.insert(*l);
				}
//				debug<<"AB\n";

				vector<IntegerVector> observedExponents;
				for(set<IntegerVector>::const_iterator a=observedExponentsSet.begin();a!=observedExponentsSet.end();a++)observedExponents.push_back(*a);

//				debug<<"AC\n";
				int numberOfObservedExponents=observedExponents.size();
				int width=(M+1)*numberOfObservedExponents+numberOfSamples;
				int height=numberOfObservedExponents*numberOfSamples;

//				debug<<"AD\n";
				FieldMatrix A(r.getField(),height,width);
//				debug<<A<<"\n";
				int K=0;
				for(vector<pair<FieldElement,Polynomial> >::const_iterator k=samples.begin();k!=samples.end();k++,K++)
				{
					for(int a=0;a<numberOfObservedExponents;a++)
					{
						FieldElement power=r.getField().zHomomorphism(1);
						for(int b=0;b<M+1;b++)
							{
								A[K*numberOfObservedExponents+a][(M+1)*a+b]=power;
								power*=k->first;
							}

						A[K*numberOfObservedExponents+a][(M+1)*numberOfObservedExponents+K]=k->second.coefficientOf(observedExponents[a]);
					}
				}
//				debug<<"OBSERVED EXPONENTS:";
				IntegerVectorList setTemp;for(auto c=observedExponents.begin();c!=observedExponents.end();c++)setTemp.push_back(*c);
//				debug.printVectorList(setTemp);
//				debug<<"Multiplication terms: 1,...,"<<"x_"<<i<<"^"<<M<<"\n";

				// we now interpolate the result

//				debug<<A<<"\n";
//				debug<<A.reduceAndComputeRank()<<"\n";
				FieldMatrix ker=A.reduceAndComputeKernel();
				//ker.reduce();
//				debug<<ker<<"\n";


				for(int c=0;c<ker.getHeight();c++)
				{
					FieldVector coefs=ker[c];
					Polynomial lifted(r);
					for(int a=0;a<numberOfObservedExponents;a++)
						for(int b=0;b<M+1;b++)
						{
//							debug<<"REPLACE"<<observedExponents[a]<<i<<b<<observedExponents[a].withIthCoordinateInserted(i,b)<<"\n";
							lifted+=Term(coefs[(M+1)*a+b],Monomial(r,observedExponents[a].withIthCoordinateInserted(i,b)));
//							debug<<Term(coefs[(M+1)*a+b],Monomial(r,observedExponents[a].withIthCoordinateInserted(i,b)));
//							debug<<" accumulating ("<<i<<b<<") "<<lifted<<"\n";
						}
//					debug<<lifted<<"\n";

					bool isOK=true;
					for(PolynomialSet::const_iterator c=p.begin();c!=p.end();c++)
					{
						Polynomial q(r);
						if(lifted.divides2(*c,&q))//does divide
						{
							q=r.zero()-q;
//							debug<<"divides:"<<lifted<<" TIMES "<<q<<" == "<<(lifted*(q))<<"\n";
						}
						else
						{
//							debug<<"does not divide.\n";
							isOK=false;
						}
					}
					if(isOK)
						{
//							debug<<"RETURNING GCD OF"<<p<<":"<<lifted<<"\n";
//						debug<<p.getRing()<<p<<"gcd:"<<lifted<<"\n\n";
							return lifted;
						}
//					debug<<"TRIED TO FIND GCD OF:\n"<<p;
//					debug<<p.getRing();
				}
//				assert(0);
				// if the interpolated result divides, then return or try to find more factors.
			}
		}
		// not successful -
		return r.one();
	}
	}
	return r.zero();//We should never go here.
}
Polynomial findFactorForRandomSubstitution(PolynomialSet const &p, int i, FieldElement v)
{
	int n=p.numberOfVariablesInRing();
	PolynomialRing r2=p.getRing().subRing(IntegerVector::allOnes(n)-IntegerVector::standardVector(n,i),true);
	PolynomialSet p2=p.withIthVariableSubstituted(r2,i,v);

	{// Check if substution is bad.
		bool allAreZero=true;
		for(PolynomialSet::const_iterator j=p2.begin();j!=p2.end();j++)if(!j->isZero())allAreZero=false;
		if(allAreZero)return r2.zero();
	}
//	debug<<"RANDOMSUBST"<<p<<"BECOMES"<<p2<<"\n";
	assert(p2.size());
	return NonMonomialPolynomialGCDForZModP(p2);
//	return p.front().e;
}

Polynomial findFactorForRandomSubstitutionAllButI(PolynomialSet p, int i)
{
	int n=p.numberOfVariablesInRing();

	for(int j=0;j<n;j++)
		if(j!=i)
		{
			int n2=p.numberOfVariablesInRing();
			PolynomialRing r2=p.getRing().subRing(IntegerVector::allOnes(n2)-IntegerVector::standardVector(n2,j>i),true);
			FieldElement v=r2.getField().random();
			p=p.withIthVariableSubstituted(r2,j>i,v);
		}
//	debug<<p<<"\n";
	assert(p.size());
	return NonMonomialPolynomialGCDForZModP(p);
}


Polynomial NonMonomialPolynomialGCDForZ(PolynomialSet p)
{
	PolynomialSet newP(p.getRing());
	IntegerVector support;
	if(simplifyPolynomialsForGCD(p,newP,support))
		return NonMonomialPolynomialGCDForZ(newP).withExpandedVariables(support,p.getRing());

	p.mark(LexicographicTermOrder());
	p.scaleMarkedCoefficientsToOne();
	p.removeDuplicates();
	{
		p=simplifyPolysViaHomogeneitySpace(p);
		p.saturate();
		p.mark(LexicographicTermOrder());
		p.scaleMarkedCoefficientsToOne();
		p.removeDuplicates();

		p=simplifyPolysViaHomogeneitySpace(p);
		p.saturate();
		p.mark(LexicographicTermOrder());
		p.scaleMarkedCoefficientsToOne();
		p.removeDuplicates();
	}
	assert(p.size());
	if(p.size()==1)return p.front();
if(1)	{
		static int i;
		i++;
		if((i==1000))
		{
//			debug<<simplifyPolysViaHomogeneitySpace(p);
//			debug<<"NonMon on:"<<p.getRing()<<p<<"\n";

			i=0;
		}

	}
	// It turns out that the above preprocessing works so well that the fancy Chinese remainder stuff is not necessary
	// - or rather its implementation is not fast enough. Using Groebner bases is better:
	p.saturate();
	p.mark(LexicographicInvertedTermOrder());
	return naiveGCD(p);


	assert(p.size());
	{
		for(PolynomialSet::iterator i=p.begin();i!=p.end();i++)
		{
			bool hasIntegerCoefficients;
			do
			{
				hasIntegerCoefficients=true;
				FieldElement s;
				for(TermMap::const_iterator j=i->terms.begin();j!=i->terms.end();j++)if(!j->second.isInteger()){hasIntegerCoefficients=false;s=j->second.multiplierGivingInteger();break;}
				if(!hasIntegerCoefficients)(*i)*=s;
			}while(!hasIntegerCoefficients);
		}
	}

	FieldElement gcdofcoef=p.front().getMarked().c;
	for(PolynomialSet::const_iterator i=p.begin();i!=p.end();i++)
	{
		FieldElement temp1=p.getRing().getField().zHomomorphism(0);
		FieldElement temp2=p.getRing().getField().zHomomorphism(0);
		gcdofcoef=gcd(gcdofcoef,i->getMarked().c,temp1,temp2);
	}

	// We should use Chinese remainer reconstruction in the following.
	// But since our examples are quite simple for now, we just stict to considering one prime at a time
	const int primes[]={9851,9967,37,13};
	for(int i=0;i<sizeof(primes)/sizeof(int);i++)
	{
		int prime=primes[i];
		PolynomialRing r2(FieldZModPZ(prime),p.getRing().getNumberOfVariables());
		PolynomialSet q=p.modularRepresentative(r2);
		bool good=true;
		for(PolynomialSet::const_iterator j=q.begin();j!=q.end();j++)if(j->isZero())good=false;
		if(!good){debug<<"FAIL ON:"<<p<<"\n";continue;}
		assert(q.size());
		Polynomial theModularGCD=NonMonomialPolynomialGCDForZModP(q);

		theModularGCD=withLexLargestCoefficientScaledTo(theModularGCD,gcdofcoef.modularRepresentative(r2.getField()));
		//		debug<<"PRIME:"<<prime<<":"<<theModularGCD<<"\n";
		Polynomial lifted=theModularGCD.integralRepresentative(p.getRing());

		bool isFactor=true;
		for(PolynomialSet::const_iterator i=p.begin();i!=p.end();i++)
			if(!lifted.divides(*i))isFactor=false;

		if(isFactor) return lifted;
	}
	return p.getRing().one();
}


const unsigned int maxcachesize=5;
class GCDCache{
public:
	vector <Polynomial> data;
	void insert(Polynomial const &p)
	{
		if(p.numberOfTerms()>2)
		{
			if(data.size()<maxcachesize)
				data.push_back(p);
			else
			{
				int i=((unsigned)rand())%maxcachesize;
				data[i]=p;
			}
		}
	}
	Polynomial commonFactor(Polynomial const &a, Polynomial const &b)
	{
		for(int i=0;i<data.size();i++)
		{
			if(a.getRing()==data[i].getRing())//IS THIS CORRECT COMPARISON????????????????????
				if(data[i].divides(a)&&data[i].divides(b))return data[i];
		}
		return a.getRing().one();
	}
};

#define USEGCDCACHE 0

#if USEGCDCACHE
static GCDCache theCache;
#endif

Polynomial polynomialGCD(Polynomial const &a_, Polynomial const &b_)
{
	Polynomial a=a_;
	Polynomial b=b_;
	assert(a.getRing().getField().isRationals());

//	if(a.isMonomial())
	IntegerVector v=min(a.greatestCommonMonomialDivisor(),b.greatestCommonMonomialDivisor());
	Monomial vm1(a.getRing(),v);
//	Monomial vm2(a.getRing(),-v);
//	a*=vm2;
//	b*=vm2;
	a.saturate();
	b.saturate();

	Polynomial ret=a.getRing().one();
	if(!a.isMonomial()&&!b.isMonomial())
	{
#if USEGCDCACHE
	Polynomial someFactor=theCache.commonFactor(a,b);
	Polynomial A(a);
	Polynomial B(b);
	someFactor.divides(a,&A);
	someFactor.divides(b,&B);
//	debug<<"SOme"<<someFactor<<"\n"<<A<<"\n"<<B<<"\n"<<a<<"\n"<<b<<"\n";

	PolynomialSet l(A.getRing());
	l.push_back(A);
	l.push_back(B);
//	debug<<"COMPUTING GCD OF:\n";debug.printPolynomialSet(l);
	ret=NonMonomialPolynomialGCDForZ(l)*someFactor;

	theCache.insert(ret);
#else
	PolynomialSet l(a.getRing());
	l.push_back(a);
	l.push_back(b);
	ret=NonMonomialPolynomialGCDForZ(l);
#endif
	}

//	debug<<"GCD COMPUTED FOR :\n";
//	debug.printPolynomialSet(l);
//	debug<<ret<<"\n\n";
	return ret*vm1;
///	return NonMonomialPolynomialGCDForZ(l);
}
#endif



#if 0


Polynomial PolynomialGCDForZModP(PolynomialRing const &r, PolynomialSet p)
{
	assert(p.size());
	IntegerVector monomialFactor=p.begin()->greatestCommonMonomialDivisor();
	for(auto i=p.begin();i!=p.end();i++)monomialFactor=min(monomialFactor,i->greatestCommonMonomialDivisor());
	p.saturate();
//	for(PolynomialSet::iterator j=p.begin();j!=p.end();j++)
//  	  *j*=Monomial(r,-monomialFactor);

    return NonMonomialPolynomialGCDForZModP(r,p)*Monomial(r,monomialFactor);
}
#endif
