/*
 * rational.h
 *
 *  Created on: Dec 9, 2013
 *      Author: anders
 */

#ifndef RATIONAL_H_
#define RATIONAL_H_

#include "integer.h"

class Rational
{
	Integer num,den;
public:
	void removeFactors()
	{
		assert(!den.isZero());
		Integer f=gcd(num,den);
		if(den.sign()<0)f=-f;
		num/=f;
		den/=f;
		num.shrink();
		den.shrink();
	}

	Rational():
		num(int64(0)),den(int64(1))
	{
	}
	Rational(int64 v):
		num(v),den(1)
	{
	}

	bool isZero()const{
		return num.isZero();
	}

	friend std::ostream &operator<<(std::ostream &f, Rational const &a)
	{
		f<<a.num<<"/"<<a.den;
		return f;
	}

	Rational& operator+=(const Rational& a)
	{
//cerr<<num<<"/"<<den<<"+"<<a.num<<"/"<<a.den<<"\n";
		Integer temp=a.den*den;
		num=num*a.den+a.num*den;
		den=temp;
		removeFactors();
		return *this;
	}
	Rational& operator-=(const Rational& a)
	{
		Integer temp=a.den*den;
		num=num*a.den-a.num*den;
		den=temp;
		removeFactors();
		return *this;
	}
	Rational& operator*=(const Rational& a)
	{
		num*=a.num;
		den*=a.den;
		removeFactors();
		return *this;
	}
	friend Rational operator-(const Rational &b)
	{
		Rational ret;
//		cerr<<"NUM:"<<b.num;
		ret.num=-b.num;
		ret.den=b.den;
//cerr<<"NUM:"<<ret.num;
		return ret;
	}
	  Rational operator+(const Rational &a)const
	  {

	    Rational ret(*this);
	    ret+=a;
	    return ret;
	  }
	  Rational operator-(const Rational &a)const
	  {
	    Rational ret(*this);
	    ret-=a;
	    return ret;
	  }
	  Rational operator*(const Rational &a)const
	  {
	    Rational ret(*this);
	    ret*=a;
	    return ret;
	  }
	  Rational operator/(const Rational &a)const
	  {
	    Rational ret;
	    ret.num=num*a.den;
	    ret.den=den*a.num;
	    return ret;
	  }
	  int sign()const
	  {
		  return num.sign()*den.sign();
	  }
		string toString(bool writeIfOne=true, bool alwaysWriteSign=false)const
		{
			bool isOne=(*this-Rational(1)).isZero();
			bool isMinusOne=(*this-Rational(-1)).isZero();

		    if(!writeIfOne && isOne)
		      {
		    	if(alwaysWriteSign)return std::string("+");
		    	return std::string("");
		      }
		    if(!writeIfOne && isMinusOne)
		      return std::string("-");
		    string S=num.toString();
		    if(!(den-Integer(1)).isZero())
		      {
		    	S=S+string("/")+den.toString();
		      }
		    if(alwaysWriteSign && sign()!=-1)
		      return std::string("+")+S;
		    return S;
		}
};
#endif /* RATIONAL_H_ */
