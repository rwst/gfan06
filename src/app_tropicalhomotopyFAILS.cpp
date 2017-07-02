/*
 * app_tropicalhomotopy.cpp
 *
 *  Created on: Jan 16, 2015
 *      Author: anders
 */

#include <assert.h>
#include "parser.h"
#include "printer.h"
#include "gfanapplication.h"
#include "log.h"
#include "matrix.h"
#include "myassert.h"
#include "linalg.h"

#include "traverser.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>

static int nthreads;
static int steps;


static void print(vector<IntegerMatrix> const &tuple)
{
	for(int i=0;i<tuple.size();i++)
		debug<<tuple[i];
	debug<<"\n";
}


int degree(IntegerMatrix const &m)//assumes entries of m are non-negative
{
	  int ret=0;
	  for(int i=0;i<m.getWidth();i++)
		  ret=max(m.column(i).sum(),ret);
	  return ret;
}

IntegerMatrix simplex(int n, int d)
{
	  IntegerMatrix ret(n,n+1);
	  for(int i=0;i<n;i++)ret[i][i+1]=d;
	  return ret;
}


int cayleyConfigurationWidth(vector<IntegerMatrix> const &tuple)
{
	  int m2=0;
	  for(int i=0;i<tuple.size();i++)
		  m2+=tuple[i].getWidth();
	  return m2;
}


class TropicalHomotopyApplication : public GFanApplication
{
	SimpleOption optionProduceTable;
	SimpleOption optionVectorInput;
	IntegerOption optionCyclic;
	IntegerOption optionNoon;
	IntegerOption optionChandra;
	IntegerOption optionKatsura;
	IntegerOption optionGaukwa;
	IntegerOption optionEco;
	IntegerOption optionNThreads;
	IntegerOption optionSteps;
public:
  const char *helpText()
  {
    return "In the future this program will compute do tropical homotopy.\n";
  }
  TropicalHomotopyApplication():
	  optionProduceTable("--table","Produces table for paper."),
	  optionVectorInput("--vectorinput","Read in a list of point configurations instead of a polynomial ring and a list of polynomials."),
	  optionCyclic("--cyclic","Use cyclic-n example instead of reading input."),
	  optionNoon("--noon","Use Noonburg-n example instead of reading input."),
	  optionChandra("--chandra","Use Chandrasekhar-n example instead of reading input."),
	  optionKatsura("--katsura","Use Katsura-n example instead of reading input."),/* Note that Verschelde's mixed volumes for the Katsura examples do not match those produced by gfan. The configurations do not seem to be identical for all n. */
	  optionGaukwa("--gaukwa","Use Gaukwa-n example instead of reading input."),
	  optionEco("--eco","Use Eco-n example instead of reading input."),
	  optionNThreads("-j","Number of threads"),
	  optionSteps("-s","Number of steps", 150)
  {
	  optionProduceTable.hide();
	  registerOptions();
  }

  const char *name()
  {
    return "_tropicalhomotopy";
  }

 class InequalityComparisonResult{//actual comparison functions were moved to the InequalityTable
 public:
	  bool empty;
	  int configurationIndex;//used for finding best
	  int columnIndex;//used for finding best
 };

class Flags{
public:
	 static const bool computeDotProductInMatrix=true;
 };

/**
 * We identify six possibly different types needed with possibly varying precission:
 * 1) The entries of the circuits (or possibly their packed representation)
 * 2) The mixed volume contribution of a single cell. This is obtained from an entry of a circuit and therefore can be represented by the above type.
 * 3) The accumulated mixed volume. This will exceed the bound of the above type in many cases. Overflows are easily checked.
 * 4) The type that dotVector uses as a result when dotting with the target. (Also used in campareInequalities)
 * 5) The intermediate type for dotVector.
 * 6) The type used in compareRevLexInverted
 *
 *
 * Type 1 and 2 are the same.
 * Type 3 is typically different.
 *
 * To simplify our design:
 *  we assume that type 4 is the same as 1 and 2. This is reasonable, as we need some bound to make type 6 efficient.
 *  we use a special (longer) type for 5, as that allows to do overflow checks at the end, assuming some bound on the target.
 *  In 6, we observe that there is no accumulation taking place. Moreover, with the assumption that 4 and 1 are the same, we only need a type with double precission to do the comparisons here, and now overflow check will be required.
 *
 *
 * To conclude, we make two types. A single precision type for 1,2,4 and a double precision type for 3,5,6
 * We further need to make assumptions on the absolute value of the entries of the target vector and the number of entries in combination to ensure that dot product computations do not overflow.
 * Overflow checks are then only needed:
 *  when casting the return value of dotVector
 *  when doing the dotDivVector operation. But since bounds are known, in most cases checks are not needed
 *  when accumulating the mixed volume
 *
 *
 *  Suggested implementations:
 *   a pair of 32/64 bit integers with only the overflow checking listed above
 *   (a pair of 32/64 bit integers with overflow checking everywhere (for debugging only))
 *   a pair of gmp integers for exact precision.
 */

class MVType{
 public:
	class Divisor{
	public:
		Divisor(MVType const &a);// for exact division
	};
	MVType(){};
	MVType(MVType const &m);
	MVType(int val);
	class MVType &operator=(int a);
	MVType operator-()const;
	MVType &operator-=(MVType a);
	MVType &operator+=(MVType a);
	MVType &operator*=(MVType a);
	friend bool operator<=(MVType const &a, MVType const &b);
	bool isZero()const;
	bool isOne()const;
	bool isNonZero()const;
	bool isNegative()const;
	bool isPositive()const;
	friend MVType min(MVType const &a, MVType const &b);
	friend MVType negabs(MVType const &a);
	static MVType computeNegativeBound(MVType * __restrict__ Ai, int w);
	static MVType dotDivVector(MVType * __restrict__ a, MVType * __restrict__ b, MVType s, MVType t, MVType::Divisor denominatorDivisor,int c, MVType boundA, MVType boundB)__attribute__ ((always_inline));
			//assigning to first argument. The return negative of the return value is an upper bound on the absolute values of the produced entries
};


/*
 * The philosophy here is that if this class overflows, then the computation needs to be restarted. Therefore
 * all overflows must be caught.
 */
class mvtyp/*: public MVType*/{
 public:
	class Divisor{
	public:
		int v;
		int shift;
		int multiplicativeInverse;
		Divisor(mvtyp const &a)// for exact division
		{ // A good place to read about these tricks seems to be the book "Hacker's Delight" by Warren.
			v=a.v;
			shift=0;
			unsigned int t=v;
			assert(t);
			while(!(t&1)){t>>=	1;shift++;}
			unsigned int inverse=t;
			while(t*inverse!=1)inverse*=2-t*inverse;
			multiplicativeInverse=inverse;
		}
	};
	class Double{
	public:
		int64 v;
		Double():v(0){};
		Double &operator+=(Double a){v+=a.v;return *this;}
		Double &operator-=(Double a){v-=a.v;return *this;}
		mvtyp castToSingle()const;
	};
 private:
public:
	int v;
	friend mvtyp operator+(mvtyp const &a, mvtyp const &b){return mvtyp(a.v+b.v);}
	friend mvtyp operator-(mvtyp const &a, mvtyp const &b){return mvtyp(a.v-b.v);}
	friend mvtyp operator*(mvtyp const &a, mvtyp const &b){return mvtyp(a.v*b.v);}
 public:
	mvtyp(){};
//	operator int()const{return v;}
	mvtyp(mvtyp const &m):v(m.v){}
	mvtyp(int val):v(val){};
	class mvtyp &operator=(int a){v=a;return *this;}
	mvtyp operator-()const{mvtyp ret;ret.v=-v;return ret;}
	mvtyp &operator-=(mvtyp a){v-=a.v;return *this;}
	mvtyp &operator+=(mvtyp a){v+=a.v;return *this;}
	mvtyp &operator*=(mvtyp a){v*=a.v;return *this;}
	friend bool operator<=(mvtyp const &a, mvtyp const &b){return a.v<=b.v;}
//	friend bool operator>=(mvtyp const &a, mvtyp const &b){return a.v>=b.v;}
	bool isZero()const{return v==0;}
	bool isOne()const{return v==1;}
	bool isNonZero()const{return v!=0;}
	bool isNegative()const{return v<0;}
	bool isPositive()const{return v>0;}
	friend int determinantSign(mvtyp const &a, mvtyp const &c, mvtyp const &b, mvtyp const &d)//NOTICE MIXED ORDER. MUST WORK BY EXTENDING
	{
		int64 r=((int64)a.v)*((int64)c.v)-((int64)b.v)*((int64)d.v);
		if(r>0)return 1;
		if(r<0)return -1;
		return 0;
	}
	Double extend()const{Double ret;ret.v=v;return ret;}
	friend Double extendedMultiplication(mvtyp const &a, mvtyp const &b){Double ret;ret.v=((int64)a.v)*((int64)b.v);return ret;}
	friend Double extendedMultiplication(mvtyp const &a, int b){Double ret;ret.v=((int64)a.v)*((int64)b);return ret;}//to be removed?
	friend mvtyp min(mvtyp const &a, mvtyp const &b){return (a.v>=b.v)?b:a;}
	friend mvtyp negabs(mvtyp const &a){return min(a,-a);}
	friend mvtyp dotDiv(mvtyp const &s, mvtyp const &a, mvtyp const &t, mvtyp const &b, mvtyp::Divisor const &q)
	{
//		return mvtyp((((int64)s.v)*((int64)a.v)+((int64)t.v)*((int64)b.v))/q.v);
//		return mvtyp((((int)((((int64)s.v)*((int64)a.v)+((int64)t.v)*((int64)b.v))>>q.shift)))*q.multiplicativeInverse);
		return mvtyp(((((int)(((unsigned int64)(((int64)s.v)*((int64)a.v)+((int64)t.v)*((int64)b.v)))>>q.shift)))*q.multiplicativeInverse));
	}
	friend void dotDivAssign(mvtyp &s, mvtyp const &a, mvtyp const &t, mvtyp const &b, mvtyp::Divisor const &q)//as above but assigning to first argument. This helps the vectorizer.
	{
		s.v=((((int)(((unsigned int64)(((int64)s.v)*((int64)a.v)+((int64)t.v)*((int64)b.v)))>>q.shift)))*q.multiplicativeInverse);
//		s.v=((((int)((((int64)s.v)*((int64)a.v)+((int64)t.v)*((int64)b.v))>>q.shift)))*q.multiplicativeInverse);
	}
	/*
	 * This function can be implemented much faster if bounds on the input vectors are known.
	 */
	static int MIN(int a, int b)
	{
		return (a<b)?a:b;
	}
	static int MAX(int a, int b)
	{
		return (a>b)?a:b;
	}
	static mvtyp computeNegativeBound(mvtyp * __restrict__ Ai, int w)
	{
		mvtyp M=0;
		mvtyp m=0;
		for(int j=0;j<w;j++)
		{
			m.v=MIN(m.v,Ai[j].v);
			M.v=MAX(M.v,Ai[j].v);
		}
		return min(m,-M);
	}
	static mvtyp quickNoShiftBounded(mvtyp * __restrict__ a, mvtyp * __restrict__ b, mvtyp s, mvtyp t, mvtyp::Divisor denominatorDivisor,int c)
	{
		mvtyp *aa=a;
		mvtyp *bb=b;
		mvtyp max=0;
	    mvtyp min=0;
	    mvtyp ret=0;
	    for(int i=0;i<c;i++)
	      {
	    	aa[i].v=((s.v*denominatorDivisor.multiplicativeInverse)*aa[i].v+
	    			(t.v*denominatorDivisor.multiplicativeInverse)*bb[i].v);
	    	min.v=MIN(min.v,aa[i].v);
	    	max.v=MAX(max.v,aa[i].v);
	      }
	    if(-max<=min)min=-max;
	    ret=min;
	    return ret;
	}
	static mvtyp dotDivVector(mvtyp * __restrict__ a, mvtyp * __restrict__ b, mvtyp s, mvtyp t, mvtyp::Divisor denominatorDivisor,int c, mvtyp boundA, mvtyp boundB)__attribute__ ((always_inline))//assigning to first argument. The return negative of the return value is an upper bound on the absolute values of the produced entries
	{
		while(((s.v|t.v)&1)==0 && denominatorDivisor.shift>0){s.v>>=1;t.v>>=1;denominatorDivisor.shift--;denominatorDivisor.v>>=1;}
#if 0
	  mvtyp *aa=(mvtyp*)__builtin_assume_aligned(a, 16);
	  mvtyp *bb=(mvtyp*)__builtin_assume_aligned(b, 16);
#else
	  mvtyp *aa=a;
	  mvtyp *bb=b;
#endif
	  mvtyp max=0;
	  mvtyp min=0;

            //debug<<"s"<<s.v<<"t"<<t.v<<"div"<<denominatorDivisor.v<<"boundA"<<boundA.v<<"boundB"<<boundB.v<<"\n";

       unsigned long int positiveResultBoundTimesD=(negabs(t).v*((int64)boundA.v)+negabs(s).v*((int64)boundB.v));

       /*      bool boolA=positiveResultBoundTimesD<((((int64)0x4000)*denominatorDivisor.v)>>denominatorDivisor.shift);
       if(boolA)//check this carefully
       {//FAST VERSION
               for(int i=0;i<c;i++)
               {
                 aa[i].v=((((signed short)s.v)*((signed short)denominatorDivisor.multiplicativeInverse))*((signed short)aa[i].v)+
                          (((signed short)t.v)*((signed short)denominatorDivisor.multiplicativeInverse))*((signed short)bb[i].v))>>denominatorDivisor.shift;
                       min.v=MIN(min.v,aa[i].v);
                       max.v=MAX(max.v,aa[i].v);
               }
               if(-max<=min)min=-max;
       }
       else*/
       //if(positiveResultBoundTimesD&1)//check this carefully
         if(positiveResultBoundTimesD<((((int64)0x40000000)*denominatorDivisor.v)>>denominatorDivisor.shift))//check this carefully
         //    if(1)
       {//FAST VERSION
        	 // debug<<s.v<<" "<<t.v<<" "<<denominatorDivisor.v<<" "<<denominatorDivisor.shift<<"\n";
        	 if(denominatorDivisor.shift==0)return quickNoShiftBounded(a,b,s,t,denominatorDivisor,c);
        	 for(int i=0;i<c;i++)
               {
                  aa[i].v=((s.v*denominatorDivisor.multiplicativeInverse)*aa[i].v+
                                               (t.v*denominatorDivisor.multiplicativeInverse)*bb[i].v)>>denominatorDivisor.shift;
                  min.v=MIN(min.v,aa[i].v);
                       max.v=MAX(max.v,aa[i].v);
               }
               if(-max<=min)min=-max;
               return min;

       }
         else
       {
         assert(0);
               for(int i=0;i<c;i++)
               {
                       dotDivAssign(aa[i],s,t,bb[i],denominatorDivisor);
                       if(max<=aa[i])max=aa[i];
                       if(aa[i]<=min)min=aa[i];
               }
               if(-max<=min)min=-max;
               }
       return min;

#if 0 //the code below is supposed to work, but does not vectorize. The code above is bad.
               //	  debug<<"s"<<s.v<<"t"<<t.v<<"boundA"<<boundA.v<<"boundB"<<boundB.v<<"\n";

	  unsigned long int positiveResultBound=(negabs(t).v*((int64)boundA.v)+negabs(s).v*((int64)boundB.v))/denominatorDivisor.v+1;

	  if(positiveResultBound<(0x40000000>>denominatorDivisor.shift))//check this carefully
	  {//FAST VERSION
		  for(int i=0;i<c;i++)
		  {
			  aa[i]=((s.v*denominatorDivisor.multiplicativeInverse)*aa[i].v+
					  (t.v*denominatorDivisor.multiplicativeInverse)*bb[i].v)>>denominatorDivisor.shift;
			  max.v=MAX(max.v,aa[i].v);
//			  if(max.v<=aa[i].v)max=aa[i].v;
			  if(min.v>=aa[i].v)min=aa[i].v;
		  }
		  if(-max<=min)min=-max;
	  }
	  else
	  {
		  for(int i=0;i<c;i++)
		  {
			  dotDivAssign(aa[i],s,t,bb[i],denominatorDivisor);
			  if(max<=aa[i])max=aa[i];
			  if(min>=aa[i])min=aa[i];
		  }
		  if(-max<=min)min=-max;
	  }
	  return min;
#endif
	  }

};




	class SingleTropicalHomotopyTraverser{
	class InequalityTable //circuit table // This table has been moved inside the IntegersectionTraverser simply because it is used nowhere else and is specific to mixed cells in Cayley configurations.
	 {
		vector<IntegerMatrix> tuple;
		vector<int> offsets;
		vector<pair<int,int> > choices;
		Matrix<mvtyp> A;//one row for each potential inequality, to store entries with indices in chosen
		Vektor<mvtyp> tempA;
		Vektor<mvtyp> Abounds;// a negative bound for each row of A, bounding the absolute value of the rows;
		vector<mvtyp> svec;//used locally
		int subconfigurationIndex;//maybe this variable is not needed
		mvtyp denominator;
		int m;
		int k;
		void printTuple()const
		{
			for(int a=0;a<tuple.size();a++)debug<<tuple[a];
		}
		bool isLegalIndex(int subconfigurationIndex, int columnIndex)const
		{
			return choices[subconfigurationIndex].first!=columnIndex && choices[subconfigurationIndex].second!=columnIndex;
		}
		  mvtyp dotVector(int subconfigurationIndex, int columnIndex, IntegerVector const &target, int onlyK=-1)const
		  {//if onlyK!=-1 then only the onlyKth subconfiguration is considered
			  mvtyp::Double ret;
			  if(onlyK!=-1)
			  {
				  if(onlyK==subconfigurationIndex)
				  {
					  int i=subconfigurationIndex;
					  ret+=extendedMultiplication(A.UNCHECKEDACCESS(i,columnIndex+offsets[subconfigurationIndex]),target.UNCHECKEDACCESS((choices)[i].second+offsets[i]));
					  ret-=extendedMultiplication(A.UNCHECKEDACCESS(i,columnIndex+offsets[subconfigurationIndex]),target.UNCHECKEDACCESS((choices)[i].first+offsets[i]));
					  ret-=extendedMultiplication(denominator,target.UNCHECKEDACCESS((choices)[i].first+offsets[i]));// the multiplication can be merged with multiplication above except that that could cause and overflow.
					  ret+=extendedMultiplication(denominator,target.UNCHECKEDACCESS(columnIndex+offsets[i]));
					  return ret.castToSingle();
				  }
				  else
				  {
					  int i=onlyK;
					  if(target.UNCHECKEDACCESS((choices)[i].first+offsets[i]))
					  {
						  ret+=extendedMultiplication(A.UNCHECKEDACCESS(i,columnIndex+offsets[subconfigurationIndex]),target.UNCHECKEDACCESS((choices)[i].second+offsets[i]));
						  ret-=extendedMultiplication(A.UNCHECKEDACCESS(i,columnIndex+offsets[subconfigurationIndex]),target.UNCHECKEDACCESS((choices)[i].first+offsets[i]));
					  }
					  return ret.castToSingle();
				  }
			  }
			  for(int i=0;i<tuple.size();i++)
			  {
				  ret+=extendedMultiplication(A.UNCHECKEDACCESS(i,columnIndex+offsets[subconfigurationIndex]),target.UNCHECKEDACCESS((choices)[i].second+offsets[i]));
				  if(i==subconfigurationIndex)
				  {
					  ret-=extendedMultiplication(A.UNCHECKEDACCESS(i,columnIndex+offsets[subconfigurationIndex]),target.UNCHECKEDACCESS((choices)[i].first+offsets[i]));
					  ret-=extendedMultiplication(denominator,target.UNCHECKEDACCESS((choices)[i].first+offsets[i]));// the multiplication can be merged with multiplication above except that that could cause and overflow.
					  ret+=extendedMultiplication(denominator,target.UNCHECKEDACCESS(columnIndex+offsets[i]));
				  }
				  else
				  {
					  ret-=extendedMultiplication(A.UNCHECKEDACCESS(i,columnIndex+offsets[subconfigurationIndex]),target.UNCHECKEDACCESS((choices)[i].first+offsets[i]));
				  }
			  }
			  return ret.castToSingle();
		  }
		  void assignDotProducts(IntegerVector const &target, int onlyK=-1)
		  {
			  int J=0;
			  for(int i=0;i<k;i++)
				  for(int j=0;j<tuple[i].getWidth();j++,J++)
					  A[k][J]=dotVector(i,j,target,onlyK);
		  }
		  bool isReverseLexInvertedLessThanZero(int subconfigurationIndex, int columnIndex)const __attribute__ ((inline))//As in ReverseLexicographicInvertedTermOrder. Compare against zero
		  {
			  int i;
			  int index=columnIndex+offsets[subconfigurationIndex];
			  for(i=0;i<subconfigurationIndex;i++)
				  if(A.UNCHECKEDACCESS(i,index).isNonZero())
				  {
					  if(choices[i].first<choices[i].second)
						  return A.UNCHECKEDACCESS(i,index).isNegative();
					  else
						  return A.UNCHECKEDACCESS(i,index).isPositive();
				  }

				  mvtyp a=A.UNCHECKEDACCESS(i,index);
				  {
					  int firstIndex=choices[i].first;
					  int secondIndex=choices[i].second;
					  int thirdIndex=columnIndex;
					  mvtyp firstValue=-a-denominator;
					  mvtyp secondValue=a;
					  mvtyp thirdValue=denominator;

					  // Bubble sort
					  if(secondIndex<firstIndex)
					  {
						  swap(secondIndex,firstIndex);
						  swap(secondValue,firstValue);
					  }
					  if(thirdIndex<secondIndex)
					  {
						  swap(secondIndex,thirdIndex);
						  swap(secondValue,thirdValue);
					  }
					  if(secondIndex<firstIndex)
					  {
						  swap(secondIndex,firstIndex);
						  swap(secondValue,firstValue);
					  }

					  if(firstValue.isNonZero())
						  return firstValue.isPositive();
					  if(secondValue.isNonZero())
						  return secondValue.isPositive();
					  if(thirdValue.isNonZero())
						  return thirdValue.isPositive();
				  }
				  i++;
				  for(;i<k;i++)
					  if(A.UNCHECKEDACCESS(i,index).isNonZero())
					  {
						  if(choices[i].first<choices[i].second)
							  return A.UNCHECKEDACCESS(i,index).isNegative();
						  else
							  return A.UNCHECKEDACCESS(i,index).isPositive();
					  }


				  return false;
		  }
	 public:
			void clearUnused()//remove?
			{
				for(int a=0;a<k+Flags::computeDotProductInMatrix;a++)
					for(int i=0;i<k;i++)
					{
						A[a][offsets[i]+choices[i].first]=0;//123456;
						A[a][offsets[i]+choices[i].second]=0;//123456;
					}
			}
			void computeABounds()
			{
				for(int i=0;i<A.getHeight();i++)
					Abounds[i]=mvtyp::computeNegativeBound(&(A[i][0]),A.getWidth());
			}
			void checkABounds()const//remove?
			{
				for(int i=0;i<A.getHeight();i++)
				{
					mvtyp M=0;
					mvtyp m=0;
					for(int j=0;j<A.getWidth();j++)
					{
						if(M<=A[i][j])M=A[i][j];
						if(A[i][j]<=m)m=A[i][j];
					}
					assert(Abounds[i]<=m);
					assert(Abounds[i]<=-M);
				}
			}
		  mvtyp getCoordinateOfInequality(int subconfigurationIndex, int columnIndex, int i, int j)const
		  {//get (i,j)th coordinate of (subconfigurationIndex,columnIndex)th inequality
			  if(i==subconfigurationIndex)
			  {
				  if(choices[i].first==j)return -A.UNCHECKEDACCESS(i,offsets[subconfigurationIndex]+columnIndex)-denominator;
				  else if(choices[i].second==j)return A.UNCHECKEDACCESS(i,offsets[subconfigurationIndex]+columnIndex);
				  else if(j==columnIndex)return denominator;
				  else return 0;

			  }
			  else
				  if(choices[i].first==j)return -A.UNCHECKEDACCESS(i,offsets[subconfigurationIndex]+columnIndex);
				  else if(choices[i].second==j)return A.UNCHECKEDACCESS(i,offsets[subconfigurationIndex]+columnIndex);
				  else return 0;
		  }
		  int sort2uniquely(int *v, int a, int b)const//a and b different
		  {
			  v[a>b]=a;
			  v[b>a]=b;
			  return 2;
		  }
		  int sort3uniquely(int *v, int a, int b, int c)const//a and b and c different
		  {
			  v[(a>b)+int(a>c)]=a;
			  v[(b>a)+int(b>c)]=b;
			  v[(c>a)+int(c>b)]=c;
			  return 3;
		  }
		  int sort4uniquely(int *v, int a, int b, int c, int d)const// a and b different and different from c and d, but c may equal d
		  {
			  if(c!=d)
			  {
			  v[(a>b)+int(a>c)+int(a>d)]=a;
			  v[(b>a)+int(b>c)+int(b>d)]=b;
			  v[(c>a)+int(c>b)+int(c>d)]=c;
			  v[(d>a)+int(d>b)+int(d>c)]=d;
			  return 4;
			  }
			  else return sort3uniquely(v,a,b,c);
		  }
		  bool compareReverseLexicographicInverted(int i1, int j1, int i2, int j2, mvtyp s1, mvtyp s2)const//s1 and s2 are always negative
		  {
			  for(int i=0;i<k;i++)
			  {
					  if(i1!=i && i2!=i)
					  {
						  int temp=determinantSign(A.UNCHECKEDACCESS(i,offsets[i1]+j1),s1,A.UNCHECKEDACCESS(i,offsets[i2]+j2),s2);
						  if(temp)
							  if(choices[i].first<choices[i].second)
								  return temp<0;
							  else
								  return temp>0;
					  }
				  int indices[4];
				  int F=choices[i].first;
				  int S=choices[i].second;
				  int toCheck;
				  if(i1==i)
					  if(i2==i)
						  toCheck=sort4uniquely(indices,F,S,j1,j2);
					  else
						  toCheck=sort3uniquely(indices,F,S,j1);
				  else
					  if(i2==i)
						  toCheck=sort3uniquely(indices,F,S,j2);
					  else
						  toCheck=sort2uniquely(indices,F,S);

				  for(int J=0;J<toCheck;J++)
				  {
					  int j=indices[J];
					  int temp=determinantSign(getCoordinateOfInequality(i1,j1,i,j),s1,getCoordinateOfInequality(i2,j2,i,j),s2);
					  if(temp>0)
						  return true;
					  else if(temp<0)
						  return false;
				  }
			  }
			  return false;
		  }
		  mvtyp getVolume()
		  {
			  return denominator;
		  }
		void replaceFirstOrSecond(bool first, int subconfigurationIndex, int newIndex, IntegerVector const &target)__attribute__ ((inline))//updates the inequality table according to replacing first or second by newIndex in the subconfigurationIndex'th configuration
		{
	//		assert(newIndex!=choices[subconfigurationIndex].first);
	//		assert(newIndex!=choices[subconfigurationIndex].second);
			int newIndex2=newIndex;for(int i=0;i<subconfigurationIndex;i++)newIndex2+=tuple[i].getWidth();
	//		debug<<"First"<<int(first)<<"subconf."<<subconfigurationIndex<<"Id"<<newIndex<<"TID"<<newIndex2<<"\n";
	//		debug<<"Before replace"<<A<<"\n";
			int oldIndex=first?choices[subconfigurationIndex].first:choices[subconfigurationIndex].second;
			(first?choices[subconfigurationIndex].first:choices[subconfigurationIndex].second)=newIndex;

			mvtyp nextDenominator=first?A[subconfigurationIndex][newIndex2]+denominator:-A[subconfigurationIndex][newIndex2];
			mvtyp t=nextDenominator;
			mvtyp::Divisor denominatorDivisor(denominator);
			for(int c=0;c<k+Flags::computeDotProductInMatrix;c++)tempA.UNCHECKEDACCESS(c)=A.UNCHECKEDACCESS(c,newIndex2);

			for(int u=0;u<m;u++)
				svec[u]=first?-A.UNCHECKEDACCESS(subconfigurationIndex,u):A.UNCHECKEDACCESS(subconfigurationIndex,u);
			mvtyp svecBound=Abounds[subconfigurationIndex];

			if(first)
				for(int j=0;j<tuple[subconfigurationIndex].getWidth();j++)
					svec[offsets[subconfigurationIndex]+j]-=denominator;//overflows should be caught
			svecBound-=denominator;


			for(int a=0;a<k+Flags::computeDotProductInMatrix;a++)
				if(first||a!=subconfigurationIndex)
				Abounds.UNCHECKEDACCESS(a)=mvtyp::dotDivVector(&A.UNCHECKEDACCESS(a,0),&(svec[0]),t,tempA.UNCHECKEDACCESS(a),denominatorDivisor,m,Abounds[a],svecBound);
	//			for(int u=0;u<m;u++)
	//				*A.UNCHECKEDACCESSRESTRICT(a,u)=dotDiv(svec[u],tempA.UNCHECKEDACCESS(a),t,A.UNCHECKEDACCESS(a,u),denominatorDivisor);


			{
				for(int a=0;a<k+Flags::computeDotProductInMatrix;a++)
					{
						A[a][offsets[subconfigurationIndex]+oldIndex]=-tempA[a];
						Abounds[a]=min(Abounds[a],negabs(tempA[a]));
					}

				if(!first)
				{
					A[subconfigurationIndex][offsets[subconfigurationIndex]+oldIndex]=-denominator;
					Abounds[subconfigurationIndex]-=denominator;//overflows should be caught
				}
			}
			denominator=nextDenominator;
	//		debug<<"After replace"<<A<<"\n";

			// We clear these unused entries of A to ensure that these columns are not chosen when comparing inequalities
			for(int i=0;i<k+Flags::computeDotProductInMatrix;i++)
			{
				A.UNCHECKEDACCESS(i,offsets[subconfigurationIndex]+choices[subconfigurationIndex].first)=0;
				A.UNCHECKEDACCESS(i,offsets[subconfigurationIndex]+choices[subconfigurationIndex].second)=0;
			}
		}
		void replaceFirst(int subconfigurationIndex, int newIndex, IntegerVector const &target){replaceFirstOrSecond(true,subconfigurationIndex,newIndex,target);}
		void replaceSecond(int subconfigurationIndex, int newIndex, IntegerVector const &target){replaceFirstOrSecond(false,subconfigurationIndex,newIndex,target);}

		InequalityTable(vector<IntegerMatrix> const &tuple_, int subconfigurationIndex_):
			tempA(tuple_.size()+Flags::computeDotProductInMatrix),
			tuple(tuple_),
			choices(tuple_.size()),
			subconfigurationIndex(subconfigurationIndex_),
			offsets(tuple_.size())
		{
			k=tuple.size();
			m=0;
			for(int i=0;i<tuple.size();i++)m+=tuple[i].getWidth();
			svec.resize(m);
			A=Matrix<mvtyp>(k+Flags::computeDotProductInMatrix,m);
			{int offset=0;for(int i=0;i<tuple.size();i++){offsets[i]=offset;offset+=tuple[i].getWidth();}}
			Abounds=Vektor<mvtyp>(k+Flags::computeDotProductInMatrix);
		}
		void setChoicesInitially()
		{
			//THIS WILL ONLY WORK FOR THE STARTING CONFIGURATION
			//sets denominators,A and choices (these must have been initialized to the right sizes
			for(int i=0;i<k;i++)
				choices[i]=pair<int, int> (i+0,i+1);
			for(int i=0;i<m;i++)
				for(int j=0;j<k;j++)
					A[j][i]=0;
			//we follow the lemma in the article. Why does the proof of the lemma split into 3 cases and not just 2?
			int a=0;
			for(int i=0;i<k;i++)
				for(int gamma=0;gamma<tuple[i].getWidth();gamma++,a++)
					if(gamma>i+1)
						for(int ii=i;ii<gamma;ii++)
							A[ii][a]=-1;
					else
						for(int ii=gamma;ii<i;ii++)
							A[ii][a]=1;
			denominator=1;
			for(int i=0;i<k;i++)Abounds[i]=-1;
	//		assignDotProducts(target);??
		}
		void compareInequalities(InequalityComparisonResult &result, IntegerVector const &target, int onlyK=-1)
		{
			bool empty=true;
			int bestConfigurationIndex;
			int bestColumnIndex;
			mvtyp targetDotBest;

			for(int i=0;i<k;i++)
			{
	            Matrix<mvtyp>::const_RowRef Ak=const_cast<const Matrix<mvtyp>&>(A)[k];
	            int offsetsi=offsets[i];
	            int tupleiwidth=tuple[i].getWidth();

				if(onlyK!=-1)if(i!=onlyK)continue;
				for(int j=0;j<tupleiwidth;j++)
					if(Flags::computeDotProductInMatrix || isLegalIndex(i,j))//unused inequalities will have value 0. Therefore isLegalIndex(i,j) is not required if values are stored.
					{
						mvtyp ineqDotTarget=Flags::computeDotProductInMatrix?Ak.UNCHECKEDACCESS(offsetsi+j):dotVector(i,j,target,onlyK);
						if(ineqDotTarget.isNegative())
						{
							if(!isReverseLexInvertedLessThanZero(i,j))
							{
								if(empty||compareReverseLexicographicInverted(bestConfigurationIndex,bestColumnIndex,i,j,ineqDotTarget,targetDotBest))
								{
									targetDotBest=ineqDotTarget;
									empty=false;
									bestConfigurationIndex=i;
									bestColumnIndex=j;
								}
							}
		//					assert(!ineq.isZero());
						}
					}
			}
			result.empty=empty;
			result.configurationIndex=bestConfigurationIndex;
			result.columnIndex=bestColumnIndex;
	//		assert(denominator>0);
		}
		void setChoicesFromEarlierHomotopy(InequalityTable const &parent, mvtyp degreeScaling, IntegerVector const &target)
		{
			//sets denominators,A and choices (these must have been initialized to the right sizes
			//columns have been added to configuration this->subconfigurationIndex
			//for that reason we need to introduce new circuits. Old circuits are preserved.
			//chioices are "relative" so no update is needed.

			choices=parent.choices;
			int numberToDrop=(subconfigurationIndex!=0) ? numberToDrop=k+1 : 0;

			choices[subconfigurationIndex-1].first-=numberToDrop;
			choices[subconfigurationIndex-1].second-=numberToDrop;

			denominator=parent.denominator;
			int offsetOld=0;
			int offsetNew=0;
			for(int i=0;i<k;i++)
			{
				int localNumberToDrop=0;
				if(i==subconfigurationIndex-1)
					localNumberToDrop=numberToDrop;
				for(int a=0;a<A.getHeight()-Flags::computeDotProductInMatrix;a++)
					for(int j=0;j<parent.tuple[i].getWidth()-localNumberToDrop;j++)
						A.UNCHECKEDACCESS(a,offsetNew+j)=(a==subconfigurationIndex ? mvtyp(1) : degreeScaling)*parent.A.UNCHECKEDACCESS(a,offsetOld+j+localNumberToDrop);
				if(i==subconfigurationIndex)
				{
					mvtyp left[20];
					for(int j=parent.tuple[i].getWidth();j<tuple[i].getWidth();j++)
						left[j]=degreeScaling;
					for(int j=parent.tuple[i].getWidth();j<tuple[i].getWidth();j++)
						for(int a=0;a<A.getHeight()-Flags::computeDotProductInMatrix;a++)
							A.UNCHECKEDACCESS(a,offsetNew+j)=0;
					for(int j=parent.tuple[i].getWidth();j<tuple[i].getWidth();j++)
					{
						for(int b=0;b<k;b++)
						{
							if(choices[subconfigurationIndex].second==b+1)
								{
									mvtyp c=tuple[i].UNCHECKEDACCESS(b,j);
									A.UNCHECKEDACCESS(subconfigurationIndex,offsetNew+j)-=c*denominator;
									Abounds.UNCHECKEDACCESS(subconfigurationIndex)+=negabs(c*denominator);
								}
						}
if(1){/*						for(int a=0;a<k+Flags::computeDotProductInMatrix;a++)
							debug<<A.UNCHECKEDACCESS(a,offsetNew+j).v<<" ";
						debug<<"\n";
*/
/*						for(int a=0;a<1;a++)
						{
															mvtyp::Double tempDouble=A.UNCHECKEDACCESS(a,offsetNew+j).extend();
								//							A.UNCHECKEDACCESS(a,offsetNew+j)=tempDouble.castToSingle();
														}
*/	/*					for(int a=0;a<k+Flags::computeDotProductInMatrix;a++)
							debug<<A.UNCHECKEDACCESS(a,offsetNew+j).v<<" ";*/
						debug<<"\n\n";
}
						for(int a=0;a<k+Flags::computeDotProductInMatrix;a++)
						{
							for(int b=0;b<k;b++)
							{
								if(choices[subconfigurationIndex].second!=b+1 &&choices[subconfigurationIndex].first!=b+1)
									{
										mvtyp c=tuple[i].UNCHECKEDACCESS(b,j);
										{
//tempDouble+=
											A.UNCHECKEDACCESS(a,offsetNew+j)+=c*A.UNCHECKEDACCESS(a,offsetNew+b+1);
											Abounds.UNCHECKEDACCESS(a)+=negabs(c*A.UNCHECKEDACCESS(a,offsetNew+b+1));
										}
									}
							}
//							A.UNCHECKEDACCESS(a,offsetNew+j)=tempDouble.castToSingle();
						}
						for(int b=0;b<k;b++)
							left[j]-=tuple[i].UNCHECKEDACCESS(b,j);
					}
					for(int j=parent.tuple[i].getWidth();j<tuple[i].getWidth();j++)
					{
						if(choices[subconfigurationIndex].second==0)
						{
							A.UNCHECKEDACCESS(subconfigurationIndex,offsetNew+j)-=left[j]*denominator;
							Abounds.UNCHECKEDACCESS(subconfigurationIndex)+=negabs(left[j]*denominator);
						}
						else if(choices[subconfigurationIndex].first!=0)
						{
							for(int a=0;a<k+Flags::computeDotProductInMatrix;a++)
							{
								A.UNCHECKEDACCESS(a,offsetNew+j)+=left[j]*A.UNCHECKEDACCESS(a,offsetNew);
								Abounds.UNCHECKEDACCESS(a)+=negabs(left[j]*denominator);
							}
						}
					}
					for(int j=0;j<parent.tuple[i].getWidth();j++)// <-- this loop has side effects on addExpressionOfCeb() above. Therefore it is placed after the loop above.
						for(int a=0;a<A.getHeight()-Flags::computeDotProductInMatrix;a++)
							A.UNCHECKEDACCESS(a,offsetNew+j)*=degreeScaling;
				}
				offsetOld+=parent.tuple[i].getWidth();
				offsetNew+=tuple[i].getWidth();
			}
			denominator*=degreeScaling;
			if(Flags::computeDotProductInMatrix)assignDotProducts(target,subconfigurationIndex);
			computeABounds();
			for(int a=0;a<k;a++)
			for(int i=0;i<k+Flags::computeDotProductInMatrix;i++)
			{
				A[i][offsets[a]+choices[a].first]=0;
				A[i][offsets[a]+choices[a].second]=0;
			}
		}
	};
		struct StackItem{
			int columnIndex;
			int configurationIndex;
			bool b;
			int choice;
			bool useFirstChanged,useSecondChanged;
			StackItem(int columnIndex_, int configurationIndex_, bool b_, int choice_, bool useFirstChanged_, bool useSecondChanged_):
				columnIndex(columnIndex_),
				configurationIndex(configurationIndex_),
				b(b_),
				choice(choice_),
				useFirstChanged(useFirstChanged_),
				useSecondChanged(useSecondChanged_)
			{
			}
		};
	public:
		vector<pair<int,int> > choices;
		IntegerVector target;
		 bool useFirstChanged;
		 bool useSecondChanged;
		 vector<StackItem> stack;
		 int eliminatedK;
		 int eliminatedKOffset;
		 vector<IntegerMatrix> tuple;
		 vector<int> offsets;
		 int m;
		 InequalityComparisonResult result;
		 InequalityTable inequalityTable;
		 void constructInequalityTableFromParent(InequalityTable const &parentTable, int degreeScaling)
		 {
			 inequalityTable.setChoicesFromEarlierHomotopy(parentTable, degreeScaling, target);
		 }
		 void constructInequalityTableInitially(int degreeScaling)
		 {
			 vector<IntegerMatrix> tempTuple;for(int i=0;i<tuple.size();i++)tempTuple.push_back(simplex(tuple.size(),1));
			 InequalityTable tempTable(tempTuple,-1);
			 tempTable.setChoicesInitially();
			 constructInequalityTableFromParent(tempTable,degreeScaling);
		 }

		 SingleTropicalHomotopyTraverser(vector<IntegerMatrix> const &tuple_, int m_, vector<pair<int,int> > const &choices_, IntegerVector const &target_, int eliminatedK_):
			 choices(choices_),
			 target(target_),
			 eliminatedK(eliminatedK_),
		     tuple(tuple_),
		     m(m_),
		     inequalityTable(tuple,eliminatedK),
		     offsets(tuple_.size())
		 {
			  eliminatedKOffset=0;for(int i=0;i<eliminatedK;i++)eliminatedKOffset+=tuple_[i].getWidth();
			  {int offset=0;for(int i=0;i<tuple.size();i++){offsets[i]=offset;offset+=tuple[i].getWidth();}}
		 }
		 virtual void process()
		{

		}
		  void print(vector<IntegerMatrix> const &tuple)
		  {
			  for(int i=0;i<tuple.size();i++)pout<<integerMatrixToFieldMatrix(tuple[i],Q);
		  }
		bool findOutgoingAndProcess(bool doProcess)//sets up useFirstChanged and useSecondChanged and processes if process is true
		{//returns true if we are at a leaf
//			inequalityTable.checkABounds();
			useFirstChanged=false;
			 useSecondChanged=false;
			 int onlyK=-1;
#if 1
			 if(eliminatedK!=-1)
				 if(target[choices[eliminatedK].first+eliminatedKOffset]==target[choices[eliminatedK].second+eliminatedKOffset])
					 onlyK=eliminatedK;
#endif

			 inequalityTable.compareInequalities(result,target,onlyK);

			 if(result.empty)
			 {
				 if(doProcess)process();
				 return true;
			 }

			 // reverse search rule:

			 mvtyp bestAtFirst=inequalityTable.getCoordinateOfInequality(result.configurationIndex,result.columnIndex,result.configurationIndex,choices[result.configurationIndex].first);
			 mvtyp bestAtSecond=inequalityTable.getCoordinateOfInequality(result.configurationIndex,result.columnIndex,result.configurationIndex,choices[result.configurationIndex].second);
			  	if(bestAtFirst.isNegative())
			  	{
			  		if(bestAtSecond.isNegative())
			  		{
			  			useFirstChanged=true;
			  			useSecondChanged=true;
			  		}
			  		else
			      	{
			      		if(bestAtSecond.isZero())
				  			useFirstChanged=true;
			      		else
			      			if(choices[result.configurationIndex].second<result.columnIndex)
					  			useFirstChanged=true;
			      	}
			  	}
			  	else
			  	{
			  		if(bestAtSecond.isNegative())
			      	{
			      		if(bestAtFirst.isZero())
			      			useSecondChanged=true;
			      		else
			      			if(choices[result.configurationIndex].first<result.columnIndex)
				      			useSecondChanged=true;
			      	}
			  		else
			      	{
			  			assert(0);
			      	}
			  	}
			  	return false;
		}
		void goToFirstChild()
		{
//			debug<<"First:  configuration index:"<<data.configurationIndex<<" column index:"<<data.columnIndex<<"\n";
			assert(useFirstChanged);
			{
				stack.push_back(StackItem(
						result.columnIndex,
						result.configurationIndex,
						false,
						choices[result.configurationIndex].first,
						useFirstChanged,
						useSecondChanged));
				choices[result.configurationIndex].first=result.columnIndex;
				inequalityTable.replaceFirst(result.configurationIndex,result.columnIndex,target);
			}
		}
		void goToSecondChild()
		{
//			debug<<"Second: configuration index:"<<data.configurationIndex<<" column index:"<<data.columnIndex<<"\n";
			assert(useSecondChanged);
			{
				stack.emplace_back(StackItem(
						result.columnIndex,
						result.configurationIndex,
						true,
						choices[result.configurationIndex].second,
						useFirstChanged,
						useSecondChanged));
				choices[result.configurationIndex].second=result.columnIndex;
				inequalityTable.replaceSecond(result.configurationIndex,result.columnIndex,target);
			}
		}
		int numberOfChildren()
		{
			return int(useFirstChanged)+int(useSecondChanged);
		}
		void goToNthChild(int n)
		{
			if(n==0)
				if(useFirstChanged)
					goToFirstChild();
				else
					goToSecondChild();
			else
				goToSecondChild();
		}
		void goBack()
		{
			StackItem &B=stack.back();
			result.columnIndex=B.columnIndex;
			result.configurationIndex=B.configurationIndex;
			if(B.b)
			{
				choices[result.configurationIndex].second=B.choice;
				inequalityTable.replaceSecond(result.configurationIndex,B.choice,target);
			}
			else
			{
				choices[result.configurationIndex].first=B.choice;
				inequalityTable.replaceFirst(result.configurationIndex,B.choice,target);
			}
			useFirstChanged=B.useFirstChanged;
			useSecondChanged=B.useSecondChanged;
			stack.pop_back();
		}
		bool atRoot()
		{
			return stack.empty();
		}
	};

  void printERR(vector<IntegerMatrix> const &tuple)
  {
	  for(int i=0;i<tuple.size();i++)debug<<integerMatrixToFieldMatrix(tuple[i],Q);
	  debug<<"{";
	  for(int i=0;i<tuple.size();i++)
		  {
			  if(i)debug<<",\n";
			  debug<<tuple[i].transposed().getRows();
		  }
	  debug<<"}\n";
  }




/*
 *  This class concatenates several homotopies to a single tree.
 */
	class TropicalRegenerationTraverser{
		// The following class is an attempt to separate homotopy data from the traversal logic.
		class Data{
		  public:
			  vector<IntegerVector> targets;
			  vector<IntegerMatrix> tuple;
			  vector<vector<IntegerMatrix> > tuples;
			  IntegerVector degrees;
			  bool isFiniteIndex(int level, int index)
			  {
				  return index>=tuple[0].getHeight()+1;
			  }

			  vector<IntegerMatrix> produceIthTuple(int i)
				{
				  int n=tuple[0].getHeight();
				  vector<IntegerMatrix> ret;
				  for(int j=0;j<tuple.size();j++)
				  {
					  if(j<i)ret.push_back(tuple[j]);
					  if(j==i)ret.push_back(combineLeftRight(simplex(n,degree(tuple[j])),tuple[j]));
					  if(j>i)ret.push_back(simplex(n,1));
				  }
				  return ret;
				}
			  Data(vector<IntegerMatrix> const &tuple_):tuple(tuple_)
			  {
				  int n=tuple[0].getHeight();

				  for(int i=0;i<tuple.size();i++)
					  degrees.push_back(degree(tuple[i]));

				  for(int i=0;i<tuple.size();i++)
					  tuples.push_back(produceIthTuple(i));

				  for(int i=0;i<tuple.size();i++)
				  {
					  IntegerVector targ;
					  for(int j=0;j<tuple.size();j++)
					  {
						  if(j==i)
							  targ=concatenation(targ,concatenation(IntegerVector::allOnes(n+1),IntegerVector(tuple[i].getWidth())));
						  else
							  targ=concatenation(targ,IntegerVector(tuples[i][j].getWidth()));
					  }
					  targets.push_back(targ);
				  }
			  };

			  vector<pair<int,int> > firstIntersection()
			  {
				  vector<pair<int,int> > ret;
				  for(int i=0;i<tuple.size();i++)
					  ret.push_back(pair<int,int>(i+0,i+1));
				  return ret;
			  }

			  void castToNextLevel(vector<pair<int,int> > const &choices, int i, int S, vector<pair<int,int> > &ret)
			  {
				  assert(ret.size()==choices.size());
				  for(int j=0;j<choices.size();j++)
					  ret[j]=choices[j];

				  assert(ret[i].first>=S);
				  assert(ret[i].second>=S);
				  ret[i].first-=S;
				  ret[i].second-=S;
			  }
		  };

	public:
		int depth;
		int counter;
		vector<SingleTropicalHomotopyTraverser> traversers;
		Data fullData;
		int level;
		bool deadEnd;
		bool isLevelLeaf;
		bool isSolutionVertex;
		vector<bool> isLevelLeafStack;
		TropicalRegenerationTraverser(vector<IntegerMatrix> const &tuple_):
			fullData(tuple_),counter(0),depth(0)
		{
			assert(tuple_.size());
			for(int i=0;i<tuple_.size();i++)
				traversers.push_back(SingleTropicalHomotopyTraverser(fullData.tuples[i],cayleyConfigurationWidth(fullData.tuples[i]),fullData.firstIntersection(),fullData.targets[i],i));
			traversers[0].constructInequalityTableInitially(fullData.degrees[0]);
			level=0;
		}
		virtual void process()
		{
			 static int i;
			 i++;
			 cerr<<"CELL"<<i<<"\n";
		}
		bool findOutgoingAndProcess(bool doProcess)
		{
			isSolutionVertex=false;
			deadEnd=false;
			isLevelLeaf=traversers[level].findOutgoingAndProcess(false);
			if(isLevelLeaf)
			{//leaf
				bool isFinite=fullData.isFiniteIndex(level,traversers[level].choices[level].first)&&fullData.isFiniteIndex(level,traversers[level].choices[level].second);
				deadEnd=!isFinite;
				if(isFinite && (level==fullData.tuple.size()-1))
				{
					isSolutionVertex=true;
					if(doProcess){process();}
					return true;
				}
			}
			return false;
		}
		int numberOfChildren()
		{
			if(isLevelLeaf&&(level==fullData.tuple.size()-1))return 0;
			if(!isLevelLeaf)
				return traversers[level].numberOfChildren();
			else
				return 1-deadEnd;
		}
		void goToNthChild(int n)
		{
			depth++;
			isLevelLeafStack.push_back(isLevelLeaf);
			if(!isLevelLeaf)
				traversers[level].goToNthChild(n);
			else
			{
				fullData.castToNextLevel(traversers[level].choices,level,fullData.tuples[level][level].getWidth()-fullData.tuples[level+1][level].getWidth(),traversers[level+1].choices);
				traversers[level+1].constructInequalityTableFromParent(traversers[level].inequalityTable,fullData.degrees[level+1]);
				level++;
			}
		}
		void print()
		{
		}
		void goBack()
		{
			depth--;
			counter++;
			deadEnd=false;
			if(traversers[level].atRoot())
				level--;
			else
				traversers[level].goBack();
			isLevelLeaf=isLevelLeafStack.back();
			isLevelLeafStack.pop_back();
		}
	};

	/*
	 * This class glues Bjarne Knudsen's to our MultiLevelIntersectionTraverser.
	 * This class should be written so that it works on any homotopy.
	 */
	class SpecializedRTraverser: public Traverser
	{
	public:
		TropicalRegenerationTraverser T;
		mvtyp::Double mixedVolume;
		int numberOfExpensiveSteps;
		SpecializedRTraverser(vector<IntegerMatrix> const &tuple_):
			 T(tuple_),
			 mixedVolume(),
			 numberOfExpensiveSteps(0)
		{
			numberOfExpensiveSteps++;
			T.findOutgoingAndProcess(false);
		}
		int  getEdgeCountNext( void )
		{
			return T.numberOfChildren();
		}

		int  moveToNext( int   index,
	                           bool  collect_info )
		{
			T.goToNthChild(index);
			numberOfExpensiveSteps++;
			T.findOutgoingAndProcess(false);
			return 0;
		}

	  void  moveToPrev( int  index )
	  {
		  T.goBack(); //index ignored
	  }

	  void  collectInfo( void )
	  {
		  if(T.isSolutionVertex)
			  mixedVolume+=T.traversers[T.level].inequalityTable.getVolume().extend();
	  }

	  void  printState( void )
	  {
		  T.print();
	  }
	};

  //will find the solutions to the tropical system using regeneration. The target coefficients are???
	vector<vector<pair<int, int> > > homotopy2(vector<IntegerMatrix> const &tuple)
	{
		assert(tuple.size());

		vector<SpecializedRTraverser> T1;
		int N=nthreads;
		if(N==0)N=1;// Even if we do not parallelize, we still need one traverser.
		T1.reserve(N);
		vector<Traverser*> I;
		for(int i=0;i<N;i++)T1.emplace_back(tuple);
		for(int i=0;i<N;i++)I.push_back(&(T1[i]));

		debug<<"Initialized.\n";

		if(N)
			traverse_threaded(&(I[0]),N,steps);
		else
			traverse_simple(I[0]);

		mvtyp::Double total;
		int totalSteps=0;
		for(int i=0;i<N;i++)
		{
			debug<<"#"<<(int)((SpecializedRTraverser*)(I[i]))->mixedVolume.v<<"Steps:"<<((SpecializedRTraverser*)(I[i]))->numberOfExpensiveSteps<<"\n";
			//debug<<((SpecializedMTraverser*)(I[i]))->T.counter;
			total+=((SpecializedRTraverser*)(I[i]))->mixedVolume;
			totalSteps+=((SpecializedRTraverser*)(I[i]))->numberOfExpensiveSteps;
		}
		debug<<"Total:"<<(int)total.v<<"\n";
		debug<<"Totalsteps:"<<totalSteps<<"\n";
		pout<<(int)total.v<<"\n";

		vector<vector<pair<int, int> > > empty;
		return empty;
	}
	void print(vector<vector<IntegerVector> > const &l)
	{
		for(int i=0;i<l.size();i++)
		{
			for(int j=0;j<l[i].size();j++)
				pout<<l[i][j];
			pout<<"\n";
		}
	}


    vector<IntegerMatrix> cyclic(int n)
	{
		vector<IntegerMatrix> ret;
		for(int i=1;i<n;i++)
		{
			IntegerMatrix m(n,n);
			for(int y=0;y<n;y++)
				for(int x=0;x<n;x++)
					m[y][x]=((x-y+n)%n)<i;
			ret.push_back(m);
		}

		IntegerVectorList c;
		c.push_back(IntegerVector::allOnes(n));
		c.push_back(IntegerVector(n));
		ret.push_back(rowsToIntegerMatrix(c,n).transposed());

		return ret;
	}
	vector<IntegerMatrix> noon(int n)
	{
		vector<IntegerMatrix> ret;
		for(int i=0;i<n;i++)
		{
			IntegerMatrix m(n,n+1);
			for(int y=0;y<n-1;y++)
				m[y+(i<=y)][y]=2;
			for(int x=0;x<n;x++)
				m[i][x]=1;
			ret.push_back(m);
		}
		return ret;
	}
	vector<IntegerMatrix> chandra(int n)
	{
		vector<IntegerMatrix> ret;
		for(int i=0;i<n;i++)
		{
			IntegerMatrix m(n,n+1);
			for(int y=0;y<n-1;y++)
				m[y][y+1]=1;
			for(int x=0;x<n;x++)
				m[i][x]+=1;
			ret.push_back(m);
		}
		return ret;
	}
	vector<IntegerMatrix> katsura(int n)
	{
		n++;
		vector<IntegerMatrix> ret;
		for(int i=0;i<n-1;i++)
		{
			IntegerMatrix m(n,n+1-((i+1)/2));
			for(int y=0;y<n-((i+1)/2);y++)
			{
				m[n-1-y][y]=1;
				m[(n-1-y-i)>0 ? (n-1-y-i) : -(n-1-y-i)][y]+=1;
			}
			m[i][m.getWidth()-1]=1;
			ret.push_back(m);
		}
		ret.push_back(combineLeftRight(IntegerMatrix::identity(n),IntegerMatrix(n,1)));
		return ret;
	}
	vector<IntegerMatrix> gaukwa(int n)
	{
		vector<IntegerMatrix> ret;
		for(int i=0;i<2*n;i++)
			ret.push_back(combineLeftRight(combineOnTop(IntegerMatrix::identity(n),i*IntegerMatrix::identity(n)),IntegerMatrix(2*n,1)));
		return ret;
	}
	vector<IntegerMatrix> eco(int n)
	{
		vector<IntegerMatrix> ret;
		for(int i=0;i<n-1;i++)
		{
			IntegerMatrix m(n,n-i);
			for(int y=0;y<n-(i+1);y++)
			{
				m[y+i][y]=1;
				m[n-1][y]=1;
				if(y)m[y-1][y]=1;
			}
			ret.push_back(m);
		}
		ret.push_back(combineLeftRight(combineOnTop(IntegerMatrix::identity(n-1),IntegerMatrix(1,n-1)),IntegerMatrix(n,1)));
		return ret;
	}
	class TimingResult{
		string name;
		string mixedVolume;
		map<int,long double> timings;
		bool isHLine;
		bool done;
		long long int n;
	public:
		TimingResult(string name_, int n_=0):
			name(name_),
			n(n_),
			done(0)
		{
			isHLine=name=="hline";
		}
		void time(int j)
		{
			if(!isHLine)
			{
				string jOption=j?"-j16":"";
				string s=string("/usr/bin/time -o gfantimingtemp -f \"%e\" gfan _tropicalhomotopy > gfanmvtemp ")+jOption+" --"+name;
				cerr<<"Running:\n"<<s<<"\n";
				system(s.c_str());
				double a;
				{
					ifstream f("gfantimingtemp");
					f>>a;
				}
				long long int mv;
				{
					ifstream f("gfanmvtemp");
					f>>mv;
				}
				mixedVolume=to_string(mv);

				timings[j]=a;
			}
		}
		void doTime(bool doSingleThreaded)
		{
			if(!done)
			{
				if(doSingleThreaded)time(0);
				time(1);
				done=true;
			}
		}
		string floatToString(double a)
		{
			stringstream s;
			s<<setw(2)<<a;
			return s.str();
		}
		string toLatexRow()
		{
			if(isHLine)return "\\hline";
			string name2=name;
			name2[0]-=32;
			return name2+"&"+to_string(n)+"&"+mixedVolume+"&"+floatToString(timings[0])+"&"+floatToString(timings[1])+"\\\\";
		}
	};
	void dumpTable(vector<TimingResult> &results)
	{
		for(int i=0;i<results.size();i++)results[i].doTime(i<results.size()-3);
		cout<<"\\begin{tabular}{|l|rrrr|}\n"<<"\\hline\n";
		cout<<"Problem&n&Mixed vol & 1 thread & 16 thr.\\\\\n";
		for(int i=0;i<results.size();i++)
			cout<<results[i].toLatexRow()<<"\n";
		cout<<"\\end{tabular}\n";
	}
	void procudeLaTeXTableForArticle()
	{
		vector<TimingResult> results;

		int leaveout=0;//1;

		results.push_back(TimingResult("hline"));
		for(int i=10;i<=16-leaveout;i++)
			results.push_back(TimingResult(string("cyclic")+to_string((long long int)i),i));
		results.push_back(TimingResult("hline"));
		dumpTable(results);
		for(int i=16;i<=23-leaveout;i++)
			results.push_back(TimingResult(string("noon")+to_string((long long int)i),i));
		results.push_back(TimingResult("hline"));
		dumpTable(results);
		for(int i=15;i<=23-leaveout;i++)
			results.push_back(TimingResult(string("chandra")+to_string((long long int)i),i));
		results.push_back(TimingResult("hline"));
		dumpTable(results);
		for(int i=15;i<=23-leaveout;i++)
			results.push_back(TimingResult(string("katsura")+to_string((long long int)i),i+1));
		results.push_back(TimingResult("hline"));
		dumpTable(results);
		for(int i=5;i<=10-leaveout;i++)
			results.push_back(TimingResult(string("gaukwa")+to_string((long long int)i),2*i));
		results.push_back(TimingResult("hline"));
		dumpTable(results);
		for(int i=19;i<=25-leaveout;i++)
			results.push_back(TimingResult(string("eco")+to_string((long long int)i),i));
		results.push_back(TimingResult("hline"));

		dumpTable(results);
	}
	int main()
	{
		if(optionProduceTable.getValue()){
			procudeLaTeXTableForArticle();
			return 0;
		}
		nthreads=optionNThreads.getValue();
		steps=optionSteps.getValue();
		vector<IntegerMatrix> tuple;
		if(optionCyclic.getValue())
			tuple=cyclic(optionCyclic.getValue());
		else if(optionNoon.getValue())
			tuple=noon(optionNoon.getValue());
		else if(optionChandra.getValue())
			tuple=chandra(optionChandra.getValue());
		else if(optionKatsura.getValue())
			tuple=katsura(optionKatsura.getValue());
		else if(optionGaukwa.getValue())
			tuple=gaukwa(optionGaukwa.getValue());
		else if(optionEco.getValue())
			tuple=eco(optionEco.getValue());
		else
		{
			IntegerVectorListList tuple1;
			if(optionVectorInput.getValue())
			{
				tuple1=FileParser(Stdin).parseIntegerVectorListList();
			}
			else
			{
				PolynomialSet g=FileParser(Stdin).parsePolynomialSetWithRing();

				tuple1=g.exponents();
			}
			for(IntegerVectorListList::const_iterator i=tuple1.begin();i!=tuple1.end();i++)
				tuple.push_back(rowsToIntegerMatrix(*i).transposed());
		}


		printERR(tuple);

		vector<vector<pair<int, int> > > completed=homotopy2(tuple);

		//    pout<<"Done\n"<<(int)completed.size()<<"\n";
		//	        printCellList(tuple,completed);

		return 0;
	}
};

//Must catch overflows to guarantee correctness
TropicalHomotopyApplication::mvtyp TropicalHomotopyApplication::mvtyp::Double::castToSingle()const//casts and checks precission
{
	assert(v<1000000000);
	assert(-v<1000000000);
	return mvtyp(v);
}

static TropicalHomotopyApplication theApplication;


#if 0
// just for testing compilers:
typedef TropicalHomotopyApplication::mvtyp mvtyp;
void inner(mvtyp * __restrict__ a, mvtyp * __restrict__ b, mvtyp s, mvtyp t, int q, mvtyp::Divisor denominatorDivisor)
{
  mvtyp *aa=(mvtyp*)__builtin_assume_aligned(a, 16);
  mvtyp *bb=(mvtyp*)__builtin_assume_aligned(b, 16);
  for(int i=0;i<256;i++)
    {
      //      aa[i]=dotDiv(aa[i],s,t,bb[i],denominatorDivisor);
      dotDivAssign(aa[i],s,t,bb[i],denominatorDivisor);
      //	    aa[i].v=((((int)(((((unsigned int64)( ((int64)s.v)*((int64)aa[i].v)+ ((int64)t.v)*((int64)bb[i].v) )))>>denominatorDivisor.shift))))*denominatorDivisor.multiplicativeInverse);
    }
}
#endif
