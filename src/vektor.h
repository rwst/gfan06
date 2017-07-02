#ifndef VEKTOR_H_INCLUDED
#define VEKTOR_H_INCLUDED

#include <vector>
#include <list>
#include <assert.h>
#include <algorithm>
#include <complex>
#include <stdio.h>

using namespace std;

typedef signed long int int64;

void outOfRange(int i, int n);

/*
 * TODO: in the process of making gfan a library to be used by POLYMAKE etc.
 * we need to:
 * - make an integer class with overflow checking. At every addition the difference
 *   between bit 31 and 32 is or'ed to a register, which can be checked later.
 * - group the current asserts into categories. Some will assert in the usual way.
 *   Others will throw exceptions. In the library version all throw exceptions.
 */

template <class typ> class Vektor{
public:
  vector<typ> v;
  int support;
public:
  //--------------
  // Constructors
  //--------------
  Vektor(const Vektor &a):v(a.v),support(a.support){
  }
  Vektor(int n):v(n){
    assert(n>=0);
    support=0;
    for(int i=0;i<n;i++)v[i]=typ(0.0);
  };
  Vektor(const typ *p, int n):v(n){
    assert(n>=0);
    support=0;for(int i=0;i<n;i++)v[i]=p[i];
  };
  ~Vektor(){
  };
  Vektor(){
    support=0;
  };
  static Vektor standardVector(int n, int i)
    {
      Vektor v(n);
      v[i]=typ(1.0);
      return v;
    }
  static Vektor allOnes(int n)
    {
      Vektor v(n);
      for(int i=0;i<n;i++)
	v[i]=typ(1.0);
      return v;
    }
  static Vektor interval(int n, int a, int b)
  {
    Vektor v(n);
    for(int i=a;i<b;i++)v[i]=1;
    return v;
  }
  //-------------
  // Conversions
  //-------------
  template<class T> explicit Vektor(const Vektor<T>& c)
	  :v(c.size())
	{
	  for(int i=0;i<size();i++)v[i]=typ(c[i]);
	}

  //--------
  // Access
  //--------
  typ& operator[](int n)
    {
      if(!(n>=0 && n<v.size()))outOfRange(n,v.size());
      return (v[n]);
    }
  const typ& operator[](int n)const{assert(n>=0 && n<v.size());return (v[n]);}
  const typ& UNCHECKEDACCESS(int n)const{return (v[n]);}
  typ& UNCHECKEDACCESS(int n){return (v[n]);}

  //-------------
  // STL related
  //-------------
  unsigned int size()const{return v.size();};
  void resize(int n){v.resize(n,0);};
  void grow(int i){if(size()<i)resize(i);}
  void push_back(typ a)
  {
    v.push_back(a);
  }
  void sort()
    {
      std::sort(v.begin(),v.end());
    }
  bool nextPermutation()
    {
      return std::next_permutation(v.begin(),v.end());
    }
  bool operator<(const Vektor & b)const
    {
      if(size()<b.size())return true;
      if(size()>b.size())return false;
      for(int i=0;i<size();i++)
	{
	  if(v[i]<b[i])return true;
	  if(b[i]<v[i])return false;
	}
      return false;
    }

  //-----------------
  // Arithmetic fast
  //-----------------
  typ sum()const{typ f=0;for(int i=0;i<size();i++)f+=v[i];return f;};
  Vektor& operator+=(const Vektor& q){assert(size()==q.size());for(int i=0;i<size();i++)v[i]+=q.v[i];return *this;}
  Vektor& operator-=(const Vektor& q){assert(size()==q.size());for(int i=0;i<size();i++)v[i]-=q.v[i];return *this;}
  inline friend typ dot(const Vektor& p, const Vektor& q){assert(p.size()==q.size());typ s=0;for(int i=0;i<p.size();i++)s+=p[i]*q[i];return s;}
 // inline friend int64 dotLong(const Vektor& p, const Vektor& q){assert(p.size()==q.size());int64 s=0;for(int i=0;i<p.size();i++)s+=(int64)p[i]*(int64)q[i];return s;}
  inline friend int64 dotLong(const Vektor& p, const Vektor& q)
  {
    assert(p.v.size()==q.v.size());
  int64 s=0;
  typename vector<typ>::const_iterator j=q.v.begin();
  for(typename std::vector<typ>::const_iterator i=p.v.begin();i!=p.v.end();i++,j++)s+=((int64)*i)*((int64)*j);
  return s;}

  inline friend void dotLong4(const Vektor& p, const Vektor& q,const Vektor& a, const Vektor& b, int64 &pa, int64 &pb, int64 &qa, int64 &qb)
  {
    assert(p.v.size()==q.v.size());
    assert(a.v.size()==b.v.size());
    assert(a.v.size()==q.v.size());
    int64 PA=0;
    int64 PB=0;
    int64 QA=0;
    int64 QB=0;
    typename vector<typ>::const_iterator qi=q.v.begin();
    typename vector<typ>::const_iterator ai=a.v.begin();
    typename vector<typ>::const_iterator bi=b.v.begin();
  for(typename std::vector<typ>::const_iterator pi=p.v.begin();pi!=p.v.end();pi++,qi++,ai++,bi++)
    {
      PA+=((int64)*pi)*((int64)*ai);
      PB+=((int64)*pi)*((int64)*bi);
      QA+=((int64)*qi)*((int64)*ai);
      QB+=((int64)*qi)*((int64)*bi);
    }
  pa=PA;
  pb=PB;
  qa=QA;
  qb=QB;
  }
  /**
   * *this represents a permutation, and we let *this act on src, and assign the result to dest.
   * If dest does not have the right size, then the routine reallocates.
   */
  inline void composeAssign(const Vektor& src, Vektor &dest)const
  {
    if(dest.size()!=src.size())dest=Vektor(src.size());
    assert(size()==src.size());
    typename std::vector<typ>::const_iterator p=v.begin();
    for(typename std::vector<typ>::iterator d=dest.v.begin();d!=dest.v.end();d++,p++)
      *d=src.UNCHECKEDACCESS(*p);
  }
  inline void composeInverseAssign(const Vektor& src, Vektor &dest)const
  {
    if(dest.size()!=src.size())dest=Vektor(src.size());
    assert(size()==src.size());
    typename std::vector<typ>::const_iterator p=v.begin();
    for(typename std::vector<typ>::const_iterator s=src.v.begin();s!=src.v.end();s++,p++)
      dest[*p]=*s;
  }

  bool operator==(const Vektor & q)const{if(size()!=q.size())return false;for(int i=0;i<size();i++)if(v[i]!=q[i])return false;return true;}
  bool operator!=(const Vektor & q)const {return !(operator==(q));}
  bool isZero() const
    {
      int n=v.size();
      for(int i=0;i<n;i++)if(v[i]!=0)return 0;
      return 1;
    }
  bool isPositive() const
    {
      int n=v.size();
      for(int i=0;i<n;i++)if(v[i]<=0)return 0;
      return 1;
    }
  bool isNonNegative() const
    {
      int n=v.size();
      for(int i=0;i<n;i++)if(v[i]<0)return 0;
      return 1;
    }
  int max()const
  {
    int ret=-0x7fffffff; //not completely correct, but kind of works for 64bit
    for(int i=0;i<v.size();i++)if(ret<v[i])ret=v[i];
    return ret;
  }
  int argMax()const
  {
	  int iret=-1;
	  int ret=-0x7fffffff; //not completely correct, but kind of works for 64bit
	  for(int i=0;i<v.size();i++)if(ret<v[i]){ret=v[i];iret=i;}
	  return iret;
  }
  int min()const
  {
    int ret=0x7fffffff;
    for(int i=0;i<v.size();i++)if(ret>v[i])ret=v[i];
    return ret;
  }
  typ infinityNorm()const
  {
	  typ a=0;
	  typ b=0;
	  for(int i=0;i<v.size();i++)
		  {
			  if(a<v[i])a=v[i];
			  if(b>v[i])b=v[i];
		  }
	  if(a>-b)return a;
	  return -b;
  }
  bool infinityNormLessThan(typ k)const
  {
	  bool aOK=true;
	  bool bOK=true;

	  for(int i=0;i<v.size();i++)
		  {
			  aOK&=(v[i]<k);
			  bOK&=(v[i]>-k);
		  }
	  return aOK&&bOK;
  }
  friend bool dependent(const Vektor& p, const Vektor& q)
    {
	  /*
      typ pp=dot(p,p);
      typ qq=dot(q,q);
      typ pq=dot(p,q);
      return pq*pq==pp*qq;
*/
	  int n=p.size();
	  assert(n==q.size());
	  int i;
	  for(i=0;i<n;i++)
	  {
		  if(p.v[i])break;
	  }
	  if(i==n)return true;
	  if(q.v[i]==0)return q.isZero();
	  int64 a=p.v[i];
	  int64 b=q.v[i];
	  for(int j=0;j<n;j++)
		  if(a*q.v[j]!=b*p.v[j])return false;
	  return true;
    }

  //-----------------
  // Arithmetic slow
  //-----------------
  inline friend Vektor operator-(const Vektor& q){return -1*q;};
  inline friend Vektor operator*(typ s, const Vektor& q){Vektor p=q;for(int i=0;i<q.size();i++)p[i]*=s;return p;}
  inline friend Vektor operator/(const Vektor& q, typ s){Vektor p=q;for(int i=0;i<q.size();i++)p[i]/=s;return p;}
  inline friend Vektor operator*(const Vektor& p, const Vektor& q){assert(p.size()==q.size());Vektor p1=p;for(int i=0;i<p.size();i++)p1.v[i]*=q.v[i];return p1;}
//  inline friend Vektor operator+(const Vektor& p, const Vektor& q){assert(p.size()==q.size());Vektor p1=p;for(int i=0;i<p.size();i++)p1[i]+=q[i];return p1;}
  inline friend Vektor operator+(const Vektor& p, const Vektor& q){if(p.size()!=q.size()){fprintf(stderr,"%i %i\n",p.size(),q.size());assert(p.size()==q.size());};Vektor p1=p;for(int i=0;i<p.size();i++)p1[i]+=q[i];return p1;}
  inline friend Vektor operator-(const Vektor& p, const Vektor& q){assert(p.size()==q.size());Vektor p1=p;for(int i=0;i<p.size();i++)p1[i]-=q[i];return p1;}
  friend Vektor max(const Vektor& p, const Vektor& q){assert(p.size()==q.size());Vektor p1=p;for(int i=0;i<p.size();i++)if(p1[i]<q[i])p1[i]=q[i];return p1;}
  friend Vektor min(const Vektor& p, const Vektor& q){assert(p.size()==q.size());Vektor p1=p;for(int i=0;i<p.size();i++)if(p1[i]>q[i])p1[i]=q[i];return p1;}
  friend Vektor coordinatewiseProduct(const Vektor& p, const Vektor& q){assert(p.size()==q.size());Vektor p1(q.size());for(int i=0;i<p.size();i++)p1[i]=p[i]*q[i];return p1;}

  //------------------
  // Monomial related
  //------------------
  int divides(const Vektor& q) const
    {
      assert(size()==q.size());
      int n=v.size();
      for(int i=0;i<n;i++)
        {
          if(v[i]>0)if(q.v[i]<v[i])return 0;
        }
      return 1;
    }
  inline friend bool relativelyPrime(const Vektor& p, const Vektor& q)
    {
      assert(p.size()==q.size());
      int n=p.size();
      for(int t=0;t<n;t++)if((p[t]>0)&&(q[t]>0)) return false;
      return true;
    }
  Vektor supportVector()const
    {
      Vektor r(v.size());
      for(int i=0;i<size();i++)
	r[i]=(v[i]!=0);
      return r;
    }

  //------------------------------
  // Subvectors and concatenation
  //------------------------------
  Vektor subvector(int begin, int end)const
    {
      assert(begin>=0);
      assert(end<=size());
      assert(end>=begin);
      Vektor ret(end-begin);
      for(int i=0;i<end-begin;i++)
	ret[i]=v[begin+i];
      return ret;
    }
  Vektor subvector(list<int> const &chosenIndices)const
  {
    Vektor ret(chosenIndices.size());
    int I=0;
    for(list<int>::const_iterator i=chosenIndices.begin();i!=chosenIndices.end();i++,I++)ret[I]=v[*i];
    return ret;
  }
  Vektor subvectorSubsetBoolean(Vektor<int> const &subset)const
  {
	  Vektor ret(subset.sum());
	  int i=0;
	  for(int j=0;j<v.size();j++)if(subset[j])ret[i++]=v[j];
	  return ret;
  }
  Vektor expandedBoolean(Vektor<int> const &subset)const // Inverse of the above, except that unknown entries are set to zero.
  {
	 Vektor ret(subset.size());
	 int i=0;
	 for(int j=0;j<subset.size();j++)
		 if(subset[j])ret[j]=v[i++];
	 return ret;
  }
  friend Vektor concatenation(Vektor const &a, Vektor const &b)
  {
    Vektor ret(a.size()+b.size());
    for(int i=0;i<a.size();i++)ret[i]=a[i];
    for(int i=0;i<b.size();i++)ret[i+a.size()]=b[i];
    return ret;
  }
  Vektor expanded(int newSize, vector<int> const &positions)
  {
	  Vektor ret(newSize);
	  assert(positions.size()==size());
	  for(int i=0;i<size();i++)ret[positions[i]]=v[i];
	  return ret;
  }
  Vektor withIthCoordinateRemoved(int i)const
  {
	  Vektor ret(size()-1);
	  for(int j=0;j<i;j++)ret[j]=v[j];
	  for(int j=i+1;j<size();j++)ret[j-1]=v[j];
	  return ret;
  }
  Vektor withIthCoordinateInserted(int i, const typ &a)const
  {
	  Vektor ret(size()+1);
	  for(int j=0;j<i;j++)ret[j]=v[j];
	  ret[i]=a;
	  for(int j=i;j<size();j++)ret[j+1]=v[j];
	  return ret;
  }

  //----------------------
  // Zero/nonZero entries
  //----------------------
  int indexOfLargestNonzeroEntry()const
  {
    int ret=-1;
    for(int i=0;i<v.size();i++)
      {
	if(v[i])ret=i;
      }
    return ret;
  }
  Vektor supportIndices()const
  {
    Vektor ret(0);
    for(int i=0;i<v.size();i++)
      if(v[i]!=0)ret.push_back(i);
    return ret;
  }
  Vektor supportAsZeroOneVector()const
  {
    Vektor ret(v.size());
    for(int i=0;i<v.size();i++)ret[i]=bool(v[i]);
    return ret;
  }
  void calcsupport(void)
    {
      support=0;
      for(int i=0;i<v.size();i++)support=(support<<1)|(((v[i]>0)==true)&1);
    }
};

typedef complex<double> ComplexNumber;

typedef Vektor<ComplexNumber> ComplexVector;
typedef Vektor<double> FloatVector;
typedef Vektor<int> IntegerVector;
typedef list<IntegerVector> IntegerVectorList;
typedef list<IntegerVectorList> IntegerVectorListList;

IntegerVectorList transposeIntegerVectorList(IntegerVectorList const &l);
IntegerVectorList multiplyIntegerVectorList(IntegerVectorList const &A, IntegerVectorList const &B);
IntegerVectorList subvectorsOfIntegerVectorList(IntegerVectorList const &l, list<int> const &chosen);
int gcdOfVector(IntegerVector const &v);
void normalizedLowLevel(IntegerVector const &v, IntegerVector &dest);
IntegerVector normalized(IntegerVector const &v);
IntegerVector randomIntegerVector(int n, int range);
/**
 * Removes duplicates and reorders list.
 */
void removeDuplicates(IntegerVectorList &l);

#endif



int gcdGFAN(int r, int s);

