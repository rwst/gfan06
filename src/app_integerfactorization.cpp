/*
 * app_primefactors.cpp
 *
 *  Created on: Nov 3, 2015
 *      Author: anders
 */

#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "reversesearch.h"
#include "termorder.h"
#include "ep_standard.h"
#include "ep_xfig.h"
#include "polyhedralcone.h"
#include "gfanapplication.h"
#include "saturation.h"
#include "field_rationals.h"
#include "field_zmodpz.h"
#include "field_rationalfunctions.h"
#include "symmetry.h"
#include "linalg.h"
#include "fieldlp.h"
#include "integer.h"
#include "polynomialgcd.h"
#include "packedmonomial.h"
//#include "gfanlib_zcone.h"

//using namespace gfan;

class IntegerFactorizationApplication : public GFanApplication
{
public:
	bool includeInDefaultInstallation()
	{
		return false;
	}
	const char *helpText()
	{
		return "TEST.\n";
	}
  IntegerFactorizationApplication()
  {
    registerOptions();
  }

  const char *name()
  {
    return "_integerfactorization";
  }

  Integer squareRootCeil(Integer N)
  {
	  Integer start=N.sqrt();
	  if(start.sign()>1)start=start-Integer(1);
	  //stupid implementation:
	  for(Integer i=start;(i-N).sign()<0;i=i+1)if((i*i-N).sign()>=0)return i;
	  assert(0);
	  return N;
  }

  Integer y(Integer N, int i, Integer sqrtceil)//returns ith square mod N
  {
	  return (Integer(i)+sqrtceil)*(Integer(i)+sqrtceil)-N;//Check overflows later
	  // returns (i+(square root of smallest square greater than or equal to sqrt(N)))^2 mod N
  }
  pair<Integer,Integer> solveQuadratic(Integer N, Integer p, Integer sqrtceil)
  {//solves (x+sqrtceil)^2=N mod p;
	  for(Integer x=(int64)0;(x-p).sign()<0;x=x+1)
	  {
		  Integer a=(x+sqrtceil)*(x+sqrtceil)-N;

		  if((((a/p)*p)-a).sign()==0)return pair<Integer,Integer>(x,mod(-sqrtceil-sqrtceil-x,p));
	  }
	  assert(0);
	  return pair<Integer,Integer>(0,0);
  }
  Integer powerModP(Integer const &base, Integer const &exponent, Integer const &p)
  {
	  if(exponent.sign()==0)return Integer(1);
	  Integer half=powerModP(base,exponent/Integer(2),p);
	  if(mod(exponent,Integer(2)).isZero())
		  return mod(half*half,p);
	  return mod(half*half*base,p);
  }
  bool isQuadraticResidueModP(Integer a, Integer p)
  {
	  return (powerModP(a,(p-Integer(1))/Integer(2),p)-Integer(1)).sign()<=0;
  }
  class Erathostenes
  {
	  int n;
	  vector<int> primes;
	  vector<int> counters;
  public:
	  Erathostenes(){n=1;}
	  int nextPrime()
	  {
		  bool ok=true;
		  do
		  {
			  ok=true;
			  n++;
			  for(int i=0;i<primes.size();i++)
			  {
				  counters[i]++;
				  if(counters[i]==primes[i])
				  {
					  ok=false;
					  counters[i]=0;
				  }
			  }
		  }
		  while(!ok);
		  primes.push_back(n);
		  counters.push_back(0);
		  return n;
	  }
	  int currentPrime(){return n;}
  };
  void printVector(vector<Integer> const &v, Printer &p)
  {
	  p<<"(";
	  for(int i=0;i<v.size();i++){if(i)p<<",";p<<v[i].toString();}
	  p<<")\n";
  }

  FieldMatrix randomMatrix(Field const &f, int height, int width)
  {
	  FieldMatrix M(f,height,width);
	  for(int i=0;i<height;i++)
		  for(int j=0;j<width;j++)
			  M[i][j]=f.zHomomorphism(rand()&1);
	  return M;
  }

  void printVectorPolynomial(vector<FieldMatrix> const &coefficients)
  {
	  for(int i=0;i<coefficients.size();i++)
		  debug<<i<<":\n"<<coefficients[i]<<"\n";
  }
  vector<FieldMatrix> productWithLinear(Field const &f, int height, int width, vector<FieldMatrix> const &coefficients, pair<FieldMatrix,FieldMatrix> const &P)
  {
	  vector<FieldMatrix> ret;
	  ret.push_back(FieldMatrix(f,height,width));
	  for(int i=0;i<coefficients.size();i++)
		  ret.push_back(coefficients[i]*P.second);
	  for(int i=0;i<coefficients.size();i++)
		  ret[i]=ret[i]+(coefficients[i]*P.first);
	  return ret;
  }
  FieldMatrix ithCoefficientOfVectorPolynomialProduct(Field const &f, int height, int width, vector<FieldMatrix> const &A, vector<FieldMatrix> const &B, int i)
  {
	  FieldMatrix ret(f,height,width);
	  for(int a=0;a<A.size();a++)
		  for(int b=0;b<B.size();b++)
		  {
			  if(a+b==i)
				  ret=ret+A[a]*B[b];
		  }

	  return ret;
  }
  int degreeOfIthColumn(vector<FieldMatrix> const &A, int i)
  {
	  int ret=-1;
	  for(int j=0;j<A.size();j++)
		  if(!A[j].transposed()[i].isZero())ret=j;
	  return ret;
  }
  int degreeOfIthColumnREV(vector<FieldMatrix> const &A, int i)
  {
	  int ret=-1;
	  for(int j=0;j<A.size();j++)
		  if(!A[j].transposed()[i].isZero()){ret=j;break;}
	  return ret;
  }
  std::pair<FieldMatrix,FieldMatrix> ALGO1(FieldMatrix Xe, IntegerVector &delta)
  {
	  int m=Xe.getHeight();
	  int n=Xe.getWidth()-Xe.getHeight();

	  debug<<"DELTA"<<delta<<"\n";

	  FieldMatrix first(Xe.getField(),m+n,m+n);
	  FieldMatrix second(Xe.getField(),m+n,m+n);

	  assert(delta.size()==m+n);

	  vector<pair<int,int> > Delta;
	  for(int i=0;i<m+n;i++)Delta.push_back(pair<int,int>(delta[i],i));
	  std::sort(Delta.begin(),Delta.end());
	  for(int i=0;i<m+n;i++)
		  first[Delta[i].second][i]=Xe.getField().zHomomorphism(1);



	  for(int i=0;i<m+n;i++)delta[i]=Delta[i].first;
	  debug<<"DELTASORTED"<<delta<<"\n";
	  Xe=Xe*first;											///!!!!!!! MISSING IN PAPER

	  IntegerVector busy(m+n);

//	  debug<<"FIRST:"<<first<<"\n";
//	  debug<<"XE:"<<Xe<<"\n";


	  for(int i=0;i<m;i++)
	  {
		  int j0;
		  for(j0=0;j0<m+n;j0++)
			  if(!Xe[i][j0].isZero() && busy[j0]==0)break;
		  if(j0<m+n)                                                 //MISSING IN PAPER??
		  {
			  busy[j0]=1;
			  for(int j=j0+1;j<m+n;j++)
			  {
				  FieldElement lambda=Xe[i][j]*Xe[i][j0].inverse();
				  for(int y=0;y<Xe.getHeight();y++)Xe[y][j]=Xe[y][j]-lambda*Xe[y][j0];
				  for(int y=0;y<first.getHeight();y++)first[y][j]=first[y][j]-lambda*first[y][j0];
			  }
		  }
		  else {
			  debug<<"NO pivot.\n";
			  assert(0);
		  }
	  }
//	  debug<<"FIRSTTT:"<<first<<"\n";
//	  debug<<"BUsy"<<busy<<"\n";

	  for(int j=0;j<m+n;j++)
		  if(busy[j])
			  for(int y=0;y<first.getHeight();y++)
			  {
				  second[y][j]=first[y][j];
				  first[y][j]=Xe.getField().zHomomorphism(0);
			  }

	  return std::pair<FieldMatrix,FieldMatrix>(first,second);
  }

  /*
   * This is an implementation of Coppersmith's block Wiedemann algorithm.
   * It follows the article "Fast computation of linear generators for matrix sequences
   * and application to the block Wiedemann algorithm" by Emmanuel Thome.
   */
  FieldMatrix findVectorInKernel(FieldMatrix const &B, int m=16, int n=16, int epsilon=1, int skip=1)
  {
	  FieldMatrix ret(B.getField(),0,B.getHeight());

	  debug<<"Height:"<<B.getHeight()<<"Width:"<<B.getWidth()<<"\n";

	  assert(B.getHeight()==B.getWidth());
	  int N=B.getHeight();

restart:
	  FieldMatrix x=randomMatrix(B.getField(),N,m);
	  FieldMatrix xT=x.transposed();
#if 0
	  FieldMatrix y=randomMatrix(B.getField(),N,n);
#else
	  FieldMatrix z=randomMatrix(B.getField(),N,n);

	  for(int i=0;i<4;i++)z=randomMatrix(B.getField(),N,n);
	  FieldMatrix y=B*z;
#endif
	  int L=N/m+N/n+epsilon;

	  debug<<"L:"<<L<<"\n";

	  vector<FieldMatrix> AX;

	  FieldMatrix Bky=y;

	  debug<<"Building A.";
	  for(int k=0;k<L;k++)
	  {
		  if(k%10==0)debug<< k<<"\n";
		  AX.push_back(xT*Bky);
		  Bky=B*Bky;
	  }
	  debug<<"L:"<<L<<"\n";
	  debug<<"Done producing polynomial A.";
//	  printVectorPolynomial(AX);

	  { //The actual algorithm:
		  IntegerVector delta(m+n);
		  int t=(m+n-1)/n+skip;
		  for(int i=0;i<delta.size();i++)delta[i]=t;

		  vector<FieldMatrix> f;
		  int ftries=0;
		  FieldMatrix mColumns(B.getField(),0,0);
		  restart2:
		  {
			  debug<<"Making random f\n";
			  f=vector<FieldMatrix>();
			  for(int i=0;i<t;i++)
				  f.push_back(combineLeftRight(randomMatrix(B.getField(),n,m),FieldMatrix(B.getField(),n,n)));
			  f.push_back(combineLeftRight(FieldMatrix(B.getField(),n,m),FieldMatrix::identity(B.getField(),n)));

			  mColumns=ithCoefficientOfVectorPolynomialProduct(B.getField(),m,m+n,AX,f,t).submatrix(0,0,m,m);
		  }
//		  debug<<mColumns<<"\n";
//		  debug<<"mColumns.rank()="<<mColumns.rank()<<" mColumns.getWidth()="<<mColumns.getWidth()<<"\n";
		  if(mColumns.rank()!=mColumns.getWidth())if(ftries++>10)goto restart; else goto restart2;

		  while(1)
		  {
			  int mintminusdeltaj=t-delta.min();
//			  debug<<"Values: mint "<<mintminusdeltaj<<" t: "<<t<<"N/m"<<N/m<<"\n";
			  if(mintminusdeltaj>N/m)break;


//			  debug<<"t-delta.max "<<mintminusdeltaj<<" N/m"<<N/m<<"\n";

#if 0
			  {//debug printing
				  vector<FieldMatrix> g;
				  debug<<"Printing figure"<<" t "<<t<<" smallest t-deltaj:"<< t-delta.max()<<"\n";
				  debug<<"AX.size()="<<(int)AX.size()<<"\n";
				  debug<<"f.size()="<<(int)f.size()<<"\n";
				  for(int i=0;i<100;i++)g.push_back(ithCoefficientOfVectorPolynomialProduct(B.getField(),m,m+n,AX,f,i));
				  for(int i=0;i<m+n;i++)
				  {
					  debug<<i<<"\tdelta="<<delta[i]<<"\t"<<"degf "<<degreeOfIthColumn(f,i)<<" "<<degreeOfIthColumnREV(f,i)<<"\t";
					  assert(degreeOfIthColumn(f,i)<=delta[i]);
					  for(int j=0;j<100;j++)
						  debug<<(g[j].transposed()[i].isZero()?" ":"*");
					  debug<<"|\n";
				  }
			  }
#endif


//			  debug<<"MULTIPLICATION\n";
//			  debug<<
			  FieldMatrix e=ithCoefficientOfVectorPolynomialProduct(B.getField(),m,m+n,AX,f,t);

			  std::pair<FieldMatrix,FieldMatrix> P=ALGO1(e,delta);
			  //			  debug<<"P:\n"<<P.first<<"\n"<<P.second<<"\n";
//			  debug<<"P:\n"<<P.second<<"\n";
//			  for(int i=0;i<f.size();i++)debug<<"f["<<i<<"]=\n"<<f[i]<<"\n";
//			  debug<<"\n";

			  f=productWithLinear(B.getField(),n,m+n,f,P);

//			  debug<<"DELTABEFORE:"<<delta<<"\n";
			  for(int i=0;i<m+n;i++)if(!P.second.transposed()[i].isZero())delta[i]++;
//			  debug<<"DELTAAFTER:"<<delta<<"\n";
			  t++;
		  }


		  FieldMatrix power=z;
		  vector<FieldMatrix> powers;				//640x64
		  for(int k=0;k<=delta.max();k++)
		  {
			  powers.push_back(power);
			  power=B*power;
		  }

		  for(int j=0;j<m+n;j++)
		  {
			  FieldMatrix w(B.getField(),B.getHeight(),1);
			  for(int k=delta[j];k>=0;k--)
			  {
				  w=w+powers[delta[j]-k]*f[k].submatrix(0,j,n,j+1);  //640x64 * 64x1
			  }

	//		  debug<<"w"<<w.transposed()<<"\n";
	//		  debug<<"B*w"<<(B*w).transposed()<<"\n";

			  if(!w.transposed()[0].isZero()&& (B*w).transposed()[0].isZero())			// ALSO CHECK CONDITION AT END OF SECTION 4 IN THE PAPER
			  {
				  debug<<"FOUND KERNEL VECTOR!\n";
				  ret.appendRow(w.transposed()[0]);
			  }
		  }
	  }
	  return ret;
  }

  vector<Integer> primeFactors(Integer N)
	{
	  int numberOfTries=3000000;//3000000;//1300000
	  int numberOfPrimes=530;//530;//230
	  vector<Integer> ret;
	  vector<Integer> factorBase;//primes
	  // choose factor base

	  Erathostenes E;

	  for(int i=0;i<numberOfPrimes;)
		  if(isQuadraticResidueModP(N,Integer(E.nextPrime()))){factorBase.push_back(Integer(E.currentPrime()));i++;}

	  printVector(factorBase,debug);

	  Integer sqrtceil=squareRootCeil(N);

	  // Build up vector of y values
	  vector<Integer> v;
	  for(int i=0;i<numberOfTries;i++)v.push_back(y(N,i,sqrtceil));

//	  printVector(v,debug);
	  for(vector<Integer>::const_iterator p=factorBase.begin();p!=factorBase.end();p++)
	  {
		  pair<Integer,Integer> pxy=solveQuadratic(N,*p,sqrtceil);
		  for(int P=pxy.first.toInt();P<v.size();P+=p->toInt())v[P]/=*p;//exact division
		  if((pxy.first-pxy.second).sign())
			  for(int P=pxy.second.toInt();P<v.size();P+=p->toInt())v[P]/=*p;//exact division
	  }
	  // Is it possible to do the exact division more than once for each P?

	  int height=0;for(int i=0;i<v.size();i++)height+=((v[i]-1).sign()==0);

	  //Now we build up matrix of Z/2Z
	  FieldZModPZ F(2);
	  FieldMatrix M(F,height,factorBase.size());

	  IntegerVector selected;

	  int I=0;
	  for(int i=0;i<v.size();i++)
		  if((v[i]-1).sign()==0)
		  {
			  for(int j=0;j<factorBase.size();j++)
			  {
				  Integer temp=(y(N,i,sqrtceil)/factorBase[j]);
				  M[I][j]=F.zHomomorphism((temp*factorBase[j]-y(N,i,sqrtceil)).sign()==0);
			  }
			  I++;
			  selected.push_back(i);
		  }

	  debug<<"Height:"<<M.getHeight()<<"Width:"<<M.getWidth()<<"\n";


	  FieldMatrix kernelVectors(F,0,M.getHeight());

	  debug<<"INPUT MATRIX: "<<M.getHeight()<<"x"<<M.getWidth()<<"\n";
	  debug<<"WE ARE LOOKING FOR COKERNEL ELEMENTS.\n";

	  //Compute kernel vector
	  if(1)
	  {
		  assert(M.getWidth()<=M.getHeight());
		  FieldMatrix B=combineLeftRight(M,FieldMatrix(M.getField(),M.getHeight(),M.getHeight()-M.getWidth())).transposed();
//		  findVectorInKernel(B);
//		  findVectorInKernel(B,16,16,1,1);
		  kernelVectors=findVectorInKernel(B,8,8,1,1);
		  kernelVectors=kernelVectors.submatrix(0,0,kernelVectors.getHeight(),M.getHeight());
	  }
	  else
	  {
		  FieldMatrix cokernel=M.transposed().reduceAndComputeKernel();
		  kernelVectors=cokernel;
	  }


	  debug<<"Number of factorisation candidates:"<<kernelVectors.getHeight()<<"\n";
	  for(int a=0;a<kernelVectors.getHeight();a++)
	  {
		  Integer s(1),s2(1);
		  for(int i=0;i<kernelVectors.getWidth();i++)if(kernelVectors[a][i].isOne())s=s*y(N,selected[i],sqrtceil);
		  for(int i=0;i<kernelVectors.getWidth();i++)if(kernelVectors[a][i].isOne())s2=s2*(Integer(selected[i])+sqrtceil);

		  Integer factor1=gcd(s2-s.sqrt(),N);
		  Integer factor2=N/factor1;
		  if((factor1-1).isZero()||(factor2-1).isZero())
		  {
			  debug<<"NOT GOOD FOR FACTORISATION\n";
		  }
		  else
		  {
			  debug<<factor1.toString()<<"*"<<factor2.toString()<<"\n";
			  ret=vector<Integer>();
			  ret.push_back(factor1);ret.push_back(factor2);
			  //return ret;
		  }
	  }
	  if(ret.size())return ret;

	  ret.push_back(N);
	  //If non-trivial

	  return ret;
	}

  int main()
  {
	  {
		  FieldZModPZ F(2);
		  FieldMatrix M(F,128,128);
		  for(int i=0;i<128;i++)
			  for(int j=0;j<128;j++)
				  M[i][j]=F.zHomomorphism((j==0 ^ i==j));
		  FieldMatrix kernelVectors(F,0,M.getWidth());
		  kernelVectors=findVectorInKernel(M,64,64,1,1);
		  debug<<kernelVectors<<"\n";
		  assert(0);
	  }





	  //Integer a="71641520761751435455133616475667090434063332228247871795429";//RSA-59
//			  Integer a="1842140345223038358851257";
		  Integer a="184214045223038358851257";
	  //	  Integer a="1842140223038358851257";
//	  Integer a="15347357";
//	  Integer a="15347";





	vector<Integer> factors=primeFactors(a);

	  printVector(factors,debug);
	  return 0;
  }
};

static IntegerFactorizationApplication theApplication;
