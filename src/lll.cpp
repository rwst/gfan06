#include "lll.h"

/*double lensq(vektor_f *a)
{
  vektor_f A1=(*a* *a);
  return A1.sum();
}
*/
int down(double f)
{
  if(f>=0)return int(f);
  return int(f-1);
}

void calcmy(IntegerMatrix const &b, FloatMatrix &my, Vektor<double> &B)
{
  FloatMatrix bs=integerToFloatMatrix(b);

  for(int k=0;k<b.getHeight();k++)
    {
      for(int j=0;j<k;j++)
	{
	  if(B[j]==0)
	    my[k][j]=0;
	  else
	    my[k][j]=dot(Vektor<double>(b[k].toVector()),bs[j])/B[j];

	  bs[k]=bs[k].toVector()-my[k][j]*bs[j].toVector();
	}
      B[k]=dot(bs[k].toVector(),bs[k].toVector());
    }
}

/*void calclambda(IntegerMatrix const &b, FloatMatrix &lambda)
{
  for(int k=0;k<b.getHeight();k++)
    {
      for(int j=0;j<=k;j++)
	{
	  double u=dot(b[k].toVector(),b[j].toVector());

	  for(int i=0;i<j;i++)
	    {
	      if(i-1==-1)
		{
		  u=lambda[i][i]*u-lambda[k][i]*lambda[j][i];
		}
	      else
		{
		  if(lambda[i-1][i-1]!=0)
		    {
		      u=lambda[i][i]*u-lambda[k][i]*lambda[j][i];
		      u/=lambda[i-1][i-1];
		    }
		}
	    }
	  lambda[k][j]=u;
	}
    }
}*/

IntegerMatrix mlll(IntegerMatrix &b, IntegerMatrix *inverseTransposedM)
{
  int k=1;
  int n=b.getHeight();// number of generators

  FloatMatrix my(n,n);//n*n

  Vektor<double> B(n);
  IntegerMatrix M=IntegerMatrix::identity(n);//n*n;
  IntegerMatrix MInverseTransposed=M;

  calcmy(b,my,B);
  while(k<n)
    {
      //size reduction
      for(int l=k-1;l>=0;l--)
	{
	  //	  calclambda(b,lambda);
	  int q=down(my[k][l]+0.5);
	  if(q)
	    {
	      b[k]=b[k].toVector()-q*b[l].toVector();
	      M[k]=M[k].toVector()-q*M[l].toVector();
	      IntegerVector temp=q*MInverseTransposed[k].toVector();
	      MInverseTransposed[l]+=temp;
	      calcmy(b,my,B);
	    }
	}
      //      calcmy(b,my,B);
      //      calclambda(b,lambda);
      //test Lovasz' condition
      if(B[k]<(0.75-my[k][k-1]*my[k][k-1])*B[k-1])
	{
	  Vektor<int> temp=b[k].toVector();b[k]=b[k-1];b[k-1]=temp;
	  Vektor<int> temp1=M[k].toVector();M[k]=M[k-1];M[k-1]=temp1;
	  temp1=MInverseTransposed[k].toVector();MInverseTransposed[k]=MInverseTransposed[k-1].toVector();MInverseTransposed[k-1]=temp1;
      calcmy(b,my,B);
	  k--;
	  if(k<1)k=1;
	}
      else k++;
    }

  if(inverseTransposedM)*inverseTransposedM=MInverseTransposed;
  return M;
}


/*int rank(basis_i B)
{
  int n=B.size();
  basis_i M=mlll(B);
  int kerdim=0;
  while((kerdim<B.size())&&(B[kerdim].iszero()))kerdim++;

  return n-kerdim;
}
*/


IntegerMatrix latticeKernelOfTransposed(IntegerMatrix const &B)
{
	debug<<"THIS FUNCTION WAS KNOWN TO FAIL on the transpose of this matrix: {"
			"(1,0,0,0),"
			"(1,0,0,0),"
			"(1,0,0,0),"
			"(0,3,0,0),"
			"(0,2,1,0),"
			"(0,1,2,0),"
			"(0,0,3,0),"
			"(0,0,0,1)} (returning just two generators) and should maybe not be used. However, the lll implementation has been repaired since then.\n";
	assert(0);
  IntegerMatrix B2=B;
  IntegerMatrix M=mlll(B2);
  int kerdim=0;
  while((kerdim<B.getHeight())&&(B2[kerdim].toVector().isZero()))
    kerdim++;

  IntegerMatrix ret(kerdim,B.getHeight());

  for(int i=0;i<kerdim;i++)
    ret[i]=M[i];

  // add asserts to check that lll reduction works correctly

  return ret;
}
