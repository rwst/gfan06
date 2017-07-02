/*
 * integer.h
 *
 *  Created on: Dec 9, 2013
 *      Author: anders
 */

#ifndef INTEGER_H_
#define INTEGER_H_

#include <gmp.h>


union PODInteger
{
	mpz_t *p;
	int64 v;
};


class Integer
{
	PODInteger data;
	void setValue(int64 v_)
	{
		data.v=(v_<<1)+1;
	}
	bool fits(int64 v)
	{
//		if(v>10)return false;
//		if(v<-10)return false;
		return ((v<<1)>>1)==v;
	}
	bool fitsIn32(int64 v)
	{
		return (int64)((((int)v)<<1)>>1)==v;
	}
	bool isImmediate()const
	{
		return data.v&1;
	}
	void extendToGmp()
	{
		int64 v=data.v>>1;
		data.p=(mpz_t*)malloc(sizeof(mpz_t));
		mpz_init(*data.p);
	    mpz_set_si(*data.p,v);
	}
public:
	void shrink()
	{
		if(!isImmediate())if(mpz_fits_sint_p(*data.p)){*this=Integer(mpz_get_si(*data.p));}
	}
	Integer()
	{
		setValue(0);
	}
	Integer(const char *s)
	{
		data.p=(mpz_t*)malloc(sizeof(mpz_t));
		mpz_init(*data.p);
		mpz_set_str(*data.p,s,10);
	}
	Integer(int64 v)
	{
		if(fits(v))
			setValue(v);
		else
		{
			data.p=(mpz_t*)malloc(sizeof(mpz_t));
			mpz_init(*data.p);
		    mpz_set_si(*data.p,v);
		}
	}
	Integer(Integer const & value_)
	{
		if(value_.isImmediate())
		{
			data.v=value_.data.v;
		}
		else
		{
			data.p=(mpz_t*)malloc(sizeof(mpz_t));
			mpz_init_set(*data.p,*value_.data.p);
		}
	}

	~Integer()
	{
		if(!isImmediate())
			{
				mpz_clear(*data.p);
				free(data.p);
			}
	}

	Integer& operator=(const Integer& a)
	    {
	      const Integer *A=(const Integer*)&a;
	      if (this != A) {
	    	  if(!isImmediate())
	    	  {
	    		  mpz_clear(*data.p);
	    		  free(data.p);
	    	  }
	    	  if(!a.isImmediate())
	    	  {
	    		  data.p=(mpz_t*)malloc(sizeof(mpz_t));
	    		  mpz_init_set(*data.p, *a.data.p);
	    	  }
	    	  else
	    		  data.v=a.data.v;
	      }
	      return *this;
	    }

	  bool isZero()const{
		  if(isImmediate())
			  return data.v==1;
		  return mpz_sgn(*data.p)==0;
	  }
	  friend std::ostream &operator<<(std::ostream &f, Integer const &a)
	  {
		if(a.isImmediate())
			f<<(a.data.v>>1);
		else
		{
			void (*freefunc)(void *, size_t);
			mp_get_memory_functions(0,0,&freefunc);
			char *str=mpz_get_str(0,10,*a.data.p);
			f<<str;
			freefunc(str,strlen(str)+1);
		}
	    return f;
	  }
	  Integer& operator+=(const Integer& a)
	    {
		  if(a.isImmediate())
		  {
			  if(isImmediate())
			  {
				  int64 res=(a.data.v>>1)+(data.v>>1);
				  if(fits(res))
					  data.v=(res<<1)+1;
				  else
					  *this=(res);
			  }
			  else
			  {
				  if((a.data.v>>1)>=0)
					  mpz_add_ui(*data.p,*data.p,a.data.v>>1);
				  else
					  mpz_sub_ui(*data.p,*data.p,-(a.data.v>>1));
			  }
		  }
		  else
		  {
			  if(isImmediate())
			  {
				  int64 temp=data.v>>1;
				  *this=a;
				  if(temp>=0)
					  mpz_add_ui(*data.p,*data.p,temp);
				  else
					  mpz_sub_ui(*data.p,*data.p,-temp);
			  }
			  else
			  {
			      mpz_add(*data.p,*data.p,*a.data.p);
			  }
		  }
	      return *this;
	    }
	  Integer& operator-=(const Integer& a)
	    {
		  if(a.isImmediate())
		  {
			  if(isImmediate())
			  {
				  int64 res=(data.v>>1)-(a.data.v>>1);
//cerr<<"\nRES:"<<res;
				  if(fits(res))
					  data.v=(res<<1)+1;
				  else
					  *this=(res);
//cerr<<*this<<"\n";
			  }
			  else
			  {
				  if((a.data.v>>1)>=0)
					  mpz_sub_ui(*data.p,*data.p,a.data.v>>1);
				  else
					  mpz_add_ui(*data.p,*data.p,-(a.data.v>>1));
			  }
		  }
		  else
		  {
			  if(isImmediate())
			  {
				  int64 temp=data.v>>1;
				  *this=a;
				  if(temp>=0)
					  mpz_sub_ui(*data.p,*data.p,temp);
				  else
					  mpz_add_ui(*data.p,*data.p,-temp);
				  mpz_neg(*data.p,*data.p);
			  }
			  else
			  {
			      mpz_sub(*data.p,*data.p,*a.data.p);
			  }
		  }
	      return *this;
	    }
	  Integer& operator*=(const Integer& a)
	    {
		  if(a.isImmediate())
		  {
			  if(isImmediate())
			  {
//				  cerr<<"HERE0\n";
				  if(fitsIn32(a.data.v>>1)&&fitsIn32(data.v>>1))
				  {
					  data.v=(((a.data.v>>1)*(data.v>>1))<<1)+1;
					  return *this;
				  }

				  int64 temp=data.v>>1;
				  data.p=(mpz_t*)malloc(sizeof(mpz_t));
				  mpz_init(*data.p);
				  mpz_set_si(*data.p,temp);

				  mpz_mul_si(*data.p,*data.p,a.data.v>>1);
//cerr<<"HERE1\n";
				  //shrink?
			  }
			  else
			  {
				  mpz_mul_si(*data.p,*data.p,a.data.v>>1);
			  }
		  }
		  else
		  {
			  if(isImmediate())
			  {
				  int64 temp=data.v>>1;
				  *this=a;
				  mpz_mul_si(*data.p,*data.p,temp);
			  }
			  else
			  {
			      mpz_mul(*data.p,*data.p,*a.data.p);
			  }
		  }

	      return *this;
	    }
	  Integer& operator/=(const Integer& a)
	    {
		  if(a.isImmediate())
		  {
			  assert(a.data.v!=1);//division by zero
			  if(isImmediate())
			  {
				  data.v=(((data.v>>1)/(a.data.v>>1))<<1)+1;
			  }
			  else
			  {
				  if((a.data.v>>1)>0)
					  mpz_div_ui(*data.p,*data.p,a.data.v>>1);
				  else
				  {
					  mpz_div_ui(*data.p,*data.p,-(a.data.v>>1));
					  mpz_mul_si(*data.p,*data.p,-1);//TODO: improve
				  }
			  }
		  }
		  else
		  {
			  if(isImmediate())
			  {
				  extendToGmp();
//				  cerr<<*this;
				  //				  if(isZero())
				  assert(0);//Not implemented yet because this is unlikely to happen????
			  }
//			  else
			  {
			      mpz_div(*data.p,*data.p,*a.data.p);
			  }
		  }

	      return *this;
	    }
	  friend Integer operator-(const Integer &b)
	  {
	    Integer ret;
	    ret-=b;
	    return ret;
	  }
	  Integer operator+(const Integer &a)const
	  {
	    Integer ret(*this);
	    ret+=a;
	    return ret;
	  }
	  Integer operator-(const Integer &a)const
	  {
	    Integer ret(*this);
	    ret-=a;
	    return ret;
	  }
	  Integer operator*(const Integer &a)const
	  {
	    Integer ret(*this);
	    ret*=a;
	    return ret;
	  }
	  Integer operator/(const Integer &a)const
	  {
	    Integer ret(*this);
	    ret/=a;
	    return ret;
	  }
#if 0
	  void madd(const Integer &a,const Integer &b)
	    {
	      mpz_t temp;
	      mpz_init(temp);
	      mpz_mul(temp,a.value,b.value);
	      mpz_add(value,value,temp);
	      mpz_clear(temp);
	    }
#endif
#if 0
	  bool operator<(const Integer &a)const
	  {

	    return mpz_cmp(value,a.value)<0;
	  }
	  bool operator==(const Integer &a)const
	  {
	    return mpz_cmp(value,a.value)==0;
	  }
#endif
	  #if 0
	  bool operator!=(const Integer &a)const
	  {
	    return mpz_cmp(value,a.value)!=0;
	  }
#endif
	  int sign()const
	  {
		  if(isImmediate())
		  {
			  if(data.v<0)return -1;
			  if((data.v>>1)>0)return 1;
			  return 0;
		  }
		  return mpz_sgn(*data.p);
	  }
#if 0
	  static Integer gcd(Integer const &a, Integer const &b, Integer &s, Integer &t)
	  {
		mpz_t r;
	    mpz_init(r);
	    mpz_gcdext(r,s.value,t.value,a.value,b.value);
	    Integer ret(r);
	    mpz_clear(r);
	    return ret;
	  }
	  /**
	   * Assigns the value to z. z must have been initialized as a gmp variable.
	   */
	  void setGmp(mpz_t z)const
	  {
	    mpz_set(z,value);
	  }
	  /**
	   * Returns a value which is useful for computing hash functions.
	   */
	  signed long int hashValue()const
	  {
	    return mpz_get_si(value);
	  }
	  bool fitsInInt()const
	  {
	    mpz_t v;
	    mpz_init(v);
	    this->setGmp(v);
	    bool ret=(mpz_fits_sint_p(v)!=0);
	    mpz_clear(v);
	    return ret;
	  }
#endif
	  int64 toInt()const //repair this function
	  {
		  if(isImmediate())
			  return data.v>>1;
		  if(mpz_fits_sint_p(*data.p))
			  return mpz_get_si(*data.p);
		  assert(0);
		  return 0;
/*		  mpz_t v;
	    mpz_init(v);
	    this->setGmp(v);
	    int ret=0;
	    if(mpz_fits_sint_p(v))
	      ret=mpz_get_si(v);
	//    else
	//      ok=false;
	    mpz_clear(v);
	    return ret;*/
	  }
	  string toString()const
	  {
		  stringstream s;
		  if(isImmediate())
		  {
//cerr<<data.v<<"\n";
			  s<<(data.v>>1);
		  }
		  else
		  {
			  void (*freefunc)(void *, size_t);
			  mp_get_memory_functions(0,0,&freefunc);
			  char *str=mpz_get_str(0,10,*data.p);
			  s<<str;
			  freefunc(str,strlen(str)+1);
		  }
		  return s.str();
	  }
	  friend Integer gcd(Integer r, Integer s)
	  {
//		  if(r.isZero())return s;
//		  if(s.isZero())return r;
		  if(r.isImmediate())
		  {
	//		  if((r.data.v>>1)==1)return 1;
			  if(s.isImmediate())
			  {
	//			  if((s.data.v>>1)==1)return 1;
				  int64 R=r.data.v>>1;
				  int64 S=s.data.v>>1;

				  if(R<0)

					  R=-R;
				    if(S<0)S=-S;

				    if(S<R)
				      {
				        int64 T;
				        T=R;
				        R=S;
				        S=T;
				      }
				    while(R!=0)
				      {
				        int64 T=S%R;
				        S=R;
				        R=T;
				      }
				    return Integer(S);
			  }
			  else
			  {
				  int64 a=r.data.v>>1;
				  if(a<0)a=-a;
				  return Integer(mpz_gcd_ui(0,*s.data.p,a));
			  }
		  }
		  else
		  {
			  if(!s.isImmediate())
			  {
				  Integer t;
				  t.data.p=(mpz_t*)malloc(sizeof(mpz_t));
				  mpz_init(*t.data.p);
				  mpz_gcd(*t.data.p,*s.data.p,*r.data.p);
				  return t;
			  }
			  else
			  {
				  if((s.data.v>>1)==1)return 1;
				  int64 a=s.data.v>>1;
				  if(a<0)a=-a;
				  return Integer(mpz_gcd_ui(0,*r.data.p,a));
			  }
		  }



		  //		  cerr<<"R"<<r<<"\n";
//		  cerr<<"S"<<s<<"\n";
//		  cerr<<r.data.v<<"\n";
//		  cerr<<s.data.v<<"\n";
		  if(r.sign()<0)
			  {
//			  cerr<<"SWAP\n";
			  r=-r;}
		    if(s.sign()<0)s=-s;

//			  cerr<<s;
		    while(!r.isZero())
		      {
		    if((s-r).sign()<0)
		      {
//				  cerr<<"---";
		        Integer t;
		        t=r;
		        r=s;
		        s=t;
		      }
	//		  cerr<<"+++";
	//		  cerr<<"R"<<r<<"\n";
	//		  cerr<<"S"<<s<<"\n";
	//			  cerr<<"R"<<r<<"\n";
	//			  cerr<<"S"<<s<<"\n";
//		        int t=s%r;
		        Integer t=s-r;
		        s=r;
		        r=t;
		      }
		    assert(!s.isZero());
//		    cerr<<s<<"\n";
//		    cerr<<"-----------\n";
//assert(0);
		    return s;
	  }
	  friend Integer mod(Integer a, Integer const &b)
	  {
		  //b must be positive
		  if(a.sign()<0)
		  {
			  a+=(-a/b+1)*b;
		  }
		  return a-(b*(a/b));
	  }
	  Integer sqrt()
	  {
		  Integer ret(*this);
		  if(ret.isImmediate())ret.extendToGmp();

		  mpz_sqrt(*ret.data.p,*ret.data.p);
		  return ret;
	  }
};

typedef PODInteger ConstInteger;

#endif /* INTEGER_H_ */
