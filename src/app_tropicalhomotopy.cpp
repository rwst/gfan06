/*
 * app_tropicalhomotopy.cpp
 *
 *  Created on: Jan 16, 2015
 *      Author: anders
 */

#include <assert.h>

#include "gfanlib_circuittableint.h"

#include "parser.h"
#include "printer.h"
#include "gfanapplication.h"
#include "log.h"
#include "myassert.h"

typedef gfan::CircuitTableInt32 mvtyp;
#include "gfanlib_tropicalhomotopy.h"
#include "gfanlib_mixedvolume.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <exception>

using namespace gfan;
using namespace gfan::MixedVolumeExamples;






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
				string s=string("/usr/bin/time -o gfantimingtemp -f \"%e\" ./gfan _tropicalhomotopy > gfanmvtemp ")+jOption+" --"+name;
				cerr<<"Running:\n"<<s<<"\n";
				int err=system(s.c_str());
				assert(err==0);
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
			s.precision(1);
			s.setf(std::ios::fixed,std::ios::floatfield);
			s<</*setw(2)<<*/a;
			return s.str();
		}
		string toLatexRow()
		{
			if(isHLine)return "\\hline";
			string name2=name;
			name2[0]-=32;
			return name2+"&"+to_string(n)+"&"+mixedVolume+"&"+floatToString(timings[0])+"&"+floatToString(timings[1])+"&"+floatToString(timings[0]/timings[1])+"\\\\";
		}
	};
	void dumpTable(vector<TimingResult> &results)
	{
		for(int i=0;i<results.size();i++)results[i].doTime(i<results.size()-3);
		cout<<"\\begin{tabular}{|l|rrrrr|}\n"<<"\\hline\n";
		cout<<"Problem&n&Mixed vol & 1 thread & 16 thr. & fac.\\\\\n";
		for(int i=0;i<results.size();i++)
			cout<<results[i].toLatexRow()<<"\n";
		cout<<"\\end{tabular}\n";
	}
	void procudeLaTeXTableForArticle()
	{
		vector<TimingResult> results;

		int leaveout=0;//0;

		results.push_back(TimingResult("hline"));
		for(int i=12;i<=18-leaveout;i++)
			results.push_back(TimingResult(string("cyclic")+to_string((long long int)i),i));
		results.push_back(TimingResult("hline"));
		dumpTable(results);
		for(int i=18;i<=26-leaveout;i++)
			results.push_back(TimingResult(string("noon")+to_string((long long int)i),i));
		results.push_back(TimingResult("hline"));
		dumpTable(results);
		for(int i=18;i<=26-leaveout;i++)
			results.push_back(TimingResult(string("chandra")+to_string((long long int)i),i));
		results.push_back(TimingResult("hline"));
		dumpTable(results);
		for(int i=15;i<=26-leaveout;i++)
			results.push_back(TimingResult(string("katsura")+to_string((long long int)i),i+1));
		results.push_back(TimingResult("hline"));
		dumpTable(results);
		for(int i=5;i<=10-leaveout;i++)
			results.push_back(TimingResult(string("gaukwa")+to_string((long long int)i),2*i));
		results.push_back(TimingResult("hline"));
		dumpTable(results);
		for(int i=19;i<=28-leaveout;i++)
			results.push_back(TimingResult(string("eco")+to_string((long long int)i),i));
		results.push_back(TimingResult("hline"));

		dumpTable(results);
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
		IntegerOption optionLoadFlow;
		IntegerOption optionNThreads;
		IntegerOption optionSteps;
		SimpleOption optionIEEE14;
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
		  optionLoadFlow("--loadflow","Load flow equations of a graph being cycle with n vertices."),
		  optionIEEE14("--ieee14","Load flow equations for the IEEE 14-bus system."),
		  optionNThreads("-j","Number of threads"),
		  optionSteps("-s","Number of steps", 500)
	  {
		  optionProduceTable.hide();
		  optionLoadFlow.hide();
		  optionIEEE14.hide();
		  registerOptions();
	  }

	  const char *name()
	  {
	    return "_tropicalhomotopy";
	  }
	  virtual bool includeInDefaultInstallation(){return false;}

	  IntMatrix rowsToIntegerMatrix(IntegerVectorList const &l)
	{
		assert(l.size());
		IntMatrix ret(l.size(),l.front().size());
		IntegerVectorList::const_iterator I=l.begin();
		for(int i=0;i<ret.getHeight();i++,I++)
			for(int j=0;j<ret.getWidth();j++)
				ret[i][j]=(*I)[j];
		return ret;
	}
	  bool isEdgeIEEE14(int i, int j, int n)
	  {
		  if(i==j)return 1;
		  const char s[]="abaebcbdbecddedidgefigijinghjkkfnmmfmllfef";
		  for(int a=0;a<sizeof(s)-1;a+=2)
			  if(i==s[a]-'a'&&j==s[a+1]-'a'||j==s[a]-'a'&&i==s[a+1]-'a')return 1;
		  return 0;
	  }
	  bool isEdge(int i, int j, int n)
	  {
		  int a=(i-j+n)%n;
		  if(a==1 || a==n-1 || a==0)return 1;
		  return 0;
	  }
	  vector<IntMatrix> edgeFlow(int n, bool special)
		{
		  vector<IntMatrix> tuple;
		  if(special)assert(n==14);
		  n--;
#if 1
			for(int swap=0;swap<2;swap++)
			for(int i=0;i<n;i++)
			{
				int nEdges=0;
				for(int k=0;k<=n;k++)nEdges+=special?isEdgeIEEE14(i+1,k,n+1):isEdge(i+1,k,n+1);
				IntMatrix A(2*n,nEdges+1);
				int K=0;
				for(int k=0;k<=n;k++)
					if(special?isEdgeIEEE14(i+1,k,n+1):isEdge(i+1,k,n+1))
					{
						A[i+n-swap*n][K]=1;
						if(k)A[k-1+swap*n][K]=1;
						K++;
					}
				tuple.push_back(A);
			}
#else
			for(int i=0;i<n;i++)
			{
				IntMatrix A(2*n,4);
				IntMatrix B(2*n,4);
				if(i!=0)B[i-1+n][0]=A[i-1][0]=1;
				if(i!=n-1)B[i+1+n][1]=A[i+1][1]=1;
				A[i+n][0]=1;
				A[i+n][1]=1;
				B[i][0]=1;
				B[i][1]=1;
				A[i][3]=A[i+n][3]=B[i][3]=B[i+n][3]=1;
				tuple.push_back(A);
				tuple.push_back(B);
			}
#endif
			for(auto i=tuple.begin();i!=tuple.end();i++)std::cerr<<*i;
			return tuple;
		}
	int main()
	{
		if(optionProduceTable.getValue()){
			procudeLaTeXTableForArticle();
			return 0;
		}
		int nthreads=optionNThreads.getValue();
		int steps=optionSteps.getValue();
		vector<IntMatrix> tuple;

		if(optionIEEE14.getValue())
			tuple=edgeFlow(14,1);
		else if(optionLoadFlow.getValue())
			tuple=edgeFlow(optionLoadFlow.getValue(),0);
		else if(optionCyclic.getValue())
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

		log1 for(auto i=tuple.begin();i!=tuple.end();i++)std::cerr<<*i;

		cout<<mixedVolume(tuple,nthreads,steps)<<"\n";

		return 0;
	}
};


static TropicalHomotopyApplication theApplication;
