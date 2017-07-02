#include "tropicalcurve.h"
#include "tropical.h"
#include "tropical2.h"
#include "buchberger.h"
#include "dimension.h"
#include "saturation.h"
#include "multiplicity.h"
#include "field_rationals.h"
#include "log.h"

bool isInTropicalVariety(PolynomialSet const &I, IntegerVector const &w, bool knownToBeHomogeneous)
{
//	debug<<"ISINON"<<I<<w<<"\n";

	for(PolynomialSet::const_iterator i=I.begin();i!=I.end();i++)
		if(initialForm(*i,w).isMonomial())return false;

/*	{
			PolynomialSet temp(I.getRing());
			for(PolynomialSet::const_iterator i=I.begin();i!=I.end();i++)
				temp.push_back(initialForm(*i,w));
			if(containsMonomial(temp))return false;
	}*/


	if(knownToBeHomogeneous)
	{
		PolynomialSet I2=I;
		WeightReverseLexicographicTermOrder T(w);
		buchberger(&I2,T,true);
		return !containsMonomial(initialForms(I2,w));
	}
	PolynomialRing R2(I.getRing().getField(),I.getRing().getNumberOfVariables()+1);
	PolynomialSet I2=I.homogenization(R2);
//	debug<<R2<<I2;
	IntegerVector w2=concatenation(w,IntegerVector(1));
	WeightReverseLexicographicTermOrder T(w2);
	buchberger(&I2,T,true);
	return !containsMonomial(initialForms(I2,w2));
}

/*
 * Let's try the following conventions:
 * The ideal must have Krull dimension 1 in the Laurent polynomial ring.
 */
IntegerVectorList tropicalCurve(PolynomialSet const &I, bool earlyExit)
{
	int stat_isInTropicalVariety=0;
	int stat_finiteLiftEasy=0;
	int stat_finiteLiftHard=0;
	log2 debug<<"tropicalCurve on:"<<I.getRing()<<I<<"\n";



	int n=I.getRing().getNumberOfVariables();
	if(n==0)
	{
		assert(0);
		IntegerVectorList ret;
	//	if(!containsMonomial(I))ret.push_back(IntegerVector(0));
		return ret;
	}


	// Find index i of a variable to project away
	PolynomialRing newRing(I.getRing().getField(),n-1);
	int i=0;


	// Move "easy" checks here



	PolynomialSet J2(newRing);

//#if 0
//	for(i=0;i<n;i++)
	{
		log2 debug<<"Doing elimination\n";
		log2 debug<<I<<"\n";
		list<int> chosenVariables;
		for(int j=0;j<n;j++)if(j!=i)chosenVariables.push_back(j);

#if 0
		IntegerVectorList A;
		A.push_back(IntegerVector::standardVector(n,i));
		A.push_back(IntegerVector::allOnes(n));
		MatrixTermOrder T(A);

		PolynomialSet J=I;
		buchberger(&J,T,true);
//		pout<<newRing<<(int)chosenVariables.size();
		J2=J.polynomialRingIntersection(newRing,&chosenVariables);
		//if(krullDimension(J2)==1)
#else
		PolynomialRing newRing2=I.getRing().withVariablesAppended("H");
		PolynomialSet J=I.homogenization(newRing2);
		IntegerVectorList A;
		A.push_back(IntegerVector::standardVector(n+1,i));
		A.push_back(IntegerVector::allOnes(n+1));
		MatrixTermOrder T(A);
		buchberger(&J,T,true);
		J.changeNumberOfVariables(I.getRing());
		J2=J.polynomialRingIntersection(newRing,&chosenVariables);
#endif


		log2 debug<<"Done doing elimination\n";
//		break;
		/*
		 * At this point we need to be careful with whether we are projecting in the torus or in affine space
		 */
	}
//#endif
	// If such variable was not found, it is because the space was already one dimensional. In this case there is an easy answer....
//	if(i==n)
	{
		if(n==1)
		{
			debug<<"Base case.\n";
//			debug<<"I:"<<I;
			int kd=krullDimension(I);
			if(kd==1)
			{
			IntegerVectorList ret;
			ret.push_back(IntegerVector::standardVector(n,0));
			if(!earlyExit)ret.push_back(-IntegerVector::standardVector(n,0));
			debug<<"Returning"<<ret<<"\n";
			return ret;
			}
			else
			{
				debug<<"Returning empty list.\n";
				return IntegerVectorList();
			}
		}
//		assert(0);
		//dimension was already 0
//		IntegerVectorList ret;
//		if(!containsMonomial(I))ret.push_back(IntegerVector(0));
//		return ret;

	}

	IntegerVectorList ret;

	// We check if the two directions which are projected to the origin are contained in the tropical variety
	for(int s=-1;s<=1;s+=2)
	{
		IntegerVector v=s*IntegerVector::standardVector(n,i);
		log2 debug<<"Tropical variety membership test of ray:"<<v<<"\n";
		stat_isInTropicalVariety++;
		if(isInTropicalVariety(I,v,false))
		{
			ret.push_back(s*IntegerVector::standardVector(n,i));
			if(earlyExit)return ret;
		}
		log2 debug<<"Done membership test.\n";
	}

	IntegerVectorList lowerDimList=tropicalCurve(J2,earlyExit);
	//assert(!lowerDimList.empty());
	/*
	 * A better approach for the following is to split lowerDimList into to lists according to whether
	 * the vector has finitely many or infinitely many lifts with respect to the prevariety of known polynomials.
	 * We then do some kind of set cover heuristics on the infinitely many ones and compute elimination polynomials.
	 * We then build up the finite number of choices.
	 * We filter everything against produced polynomials.
	 * Do initial ideal computation and collect produced polynomials for filter.
	 */


	IntegerVectorList complicatedList;
	IntegerVectorList toCheck;

	for(IntegerVectorList::const_iterator k=lowerDimList.begin();k!=lowerDimList.end();k++)
	{
		log2 debug<<"Doing polyhedral computations1\n";
#if 1
		IntegerVectorList generators;generators.push_back((concatenation(IntegerVector(1),*k)));
		IntegerVectorList lin;lin.push_back((IntegerVector::standardVector(n,i)));
		PolyhedralCone C=PolyhedralCone::givenByRays(generators,lin,n);

		PolyhedralFan F(n);
		C.canonicalize();
		F.insert(C);

//		PolyhedralFan slice(n);
//		slice.insert(C);

//		debug<<"CCCCCC"<<generators<<lin<<C;
		{
		for(PolynomialSet::const_iterator i=I.begin();i!=I.end();i++)
		{
//			PolyhedralFan temp=refinement(slice,PolyhedralFan::bergmanOfPrincipalIdeal(*i));
//			F=refinement(temp,F);
			F=refinement(PolyhedralFan::bergmanOfPrincipalIdeal(*i),F);
//			if(F.getMaxDimension()==1)break;
		}
		}
#else
		debug << "OLD" << I << "\n";
		PolynomialRing r(Q,2);
		PolynomialSet newI(r);
		for(PolynomialSet::const_iterator i2=I.begin();i2!=I.end();i2++)
		{
			Polynomial p(r);
			for(TermMap::const_iterator j=i2->terms.begin();j!=i2->terms.end();j++)
				p+=Term(Q.zHomomorphism(1),Monomial(r,j->first.exponent[i]*IntegerVector::standardVector(2,0)+dot(j->first.exponent,concatenation(IntegerVector(1),*k))*IntegerVector::standardVector(2,1)));
			newI.push_back(p);
		}
		debug << "NEW" << newI << "\n";
		assert(0);

		PolyhedralFan::fullSpace(2) Fnew=PolyhedralFan::fullSpace(2);

		for(PolynomialSet::const_iterator i=newI.begin();i!=newI.end();i++)
		{
			Fnew=refinement(PolyhedralFan::normalFanOfNewtonPolytope(*i),F);
		}

		if(Fnew.conesBegin()->linealitySpace().dimension()==1)
		{
			IntegerVectorList toCheckNew=Fnew.getRays();
			for(PolyhedralFan::coneIterator i=Fnew.begin();i!=Fnew.end();i++)
				if(i->dimension()==2)
					toCheckNew.push_back(i->getRelativeInteriorPoint());
			IntegerVectorList toCheck;
			for(IntegerVector::const_iterator i=toCheckNew.begin();i!=toCheckNew.end();i++)
				toCheck.push_back();
		}
		else
		{
			assert(0);
		}
		PolyhedralFan F(n);

#endif


//		debug<<"COOONE"<<*k<<F;

		int fanDim=F.getMaxDimension();


/*		if(fanDim==0)
		{
			PolyhedralFan F(n);
			F.insert(C);


			debug<<"DEBUG";
			debug<<"lifting:"<<*k<<"\n";


			for(PolynomialSet::const_iterator i=I.begin();i!=I.end();i++)
			{
//				debug<<"COOONE"<<*i<<F.printWithI;

				debug<<"Exponents:"<<i->exponents()<<"\n";

				F.printWithIndices(&debug,FPF_default);
				PolyhedralFan t=PolyhedralFan::bergmanOfPrincipalIdeal(*i);
				t.printWithIndices(&debug,FPF_default);
				F=refinement(t,F);
				F.printWithIndices(&debug,FPF_default);
			}
		}*/

		log2 debug<<"Done with polyhedral computations\n";

		log2 debug<<"Intersection fan dim "<<fanDim<<"\n";

		assert(fanDim!=0);
		if(fanDim==1)
		{
			stat_finiteLiftEasy++;
			IntegerVectorList temp=F.getRays(fanDim);
			IntegerVectorList candidates;
			for(IntegerVectorList::const_iterator j=temp.begin();j!=temp.end();j++)
				if(!j->subvector(1,n).isZero())
					candidates.push_back(*j);

			log2 debug<<"Number of candidates "<<(int)candidates.size()<<"\n";

			if(candidates.size()==1)
			{
				ret.push_back(candidates.front());
				if(earlyExit)return ret;
			}
			else
			{
				toCheck.splice(toCheck.end(),candidates);
			}
			assert(i==0);
		}
		else
		{
			stat_finiteLiftHard++;
			complicatedList.push_back(*k);
		}
	}



	log2 debug<<"Checking complicated list:\n";
	log2 debug<<complicatedList<<"\n";

	PolynomialSet additionalPolys=I.getRing();

	while(!complicatedList.empty())
	{
		IntegerVector supportSum(n-1);
		for(IntegerVectorList::const_iterator k=complicatedList.begin();k!=complicatedList.end();k++)
			supportSum+=k->supportAsZeroOneVector();

		log2 debug<<"Support sum:"<<supportSum<<"\n";
		int j=supportSum.argMax();


		int j2=j;
		if(j2>=i)j2++;
		list<int> chosenVariables;
		chosenVariables.push_back(i);
		chosenVariables.push_back(j2);
		IntegerVectorList A;
#if 0
		A.push_back(IntegerVector::allOnes(n)-IntegerVector::standardVector(n,i)-IntegerVector::standardVector(n,j2));
		A.push_back(IntegerVector::allOnes(n));
		MatrixTermOrder T(A);
		PolynomialSet J=I;

		debug<<"Doing elimination\n";
		buchberger(&J,T,true);
		PolynomialRing newRing2(newRing.getField(),2);
		PolynomialSet P=J.polynomialRingIntersection(newRing2,&chosenVariables);
#else
		A.push_back(IntegerVector::allOnes(n+1)-IntegerVector::standardVector(n+1,i)-IntegerVector::standardVector(n+1,j2)-IntegerVector::standardVector(n+1,n));
		A.push_back(IntegerVector::allOnes(n+1));
		MatrixTermOrder T(A);

		PolynomialRing newRing3=I.getRing().withVariablesAppended("H");
		PolynomialSet J=I.homogenization(newRing3);
		log2 debug<<"Doing elimination\n";
		buchberger(&J,T,true);
		J.changeNumberOfVariables(I.getRing());
		PolynomialRing newRing2(newRing.getField(),2);
		PolynomialSet P=J.polynomialRingIntersection(newRing2,&chosenVariables);
#endif
		additionalPolys.splice(additionalPolys.end(),J);
		log2 debug<<"Done eliminating\n";
		log2 debug<<P<<"\n";
		assert(P.size()==1);

		log2 debug<<"Doing polyhedral computation\n";
		//		debug<<"PFRONT:"<<P.front()<<P.front().exponents()<<"\n";
				//		PolyhedralFan FF=PolyhedralFan::normalFanOfNewtonPolytope(P.front());
		PolyhedralFan FF=PolyhedralFan::bergmanOfPrincipalIdeal(P.front());
		IntegerVectorList F;
		if(FF.dimensionOfLinealitySpace()==1)
			{
				assert(FF.conesBegin()!=FF.conesEnd());
				IntegerVectorList l=FF.conesBegin()->generatorsOfLinealitySpace();
				assert(l.size()==1);
				F.push_back(l.front());
				F.push_back(-l.front());
//				debug<<"LINgen"<<F<<"\n";
			}
		else
		{
//			debug<<"Calling getrays\n";
			F=FF.getRays();
		}

		log2 debug<<"Potential lifts"<<F<<"\n";


		for(IntegerVectorList::iterator k=complicatedList.begin();k!=complicatedList.end();)
			if((*k)[j]!=0)
		{
	//		debug<<I.getRing()<<"\n"<<I;
			log2 debug<<"Checking:"<<*k<<" with chosen coordinate:"<<j<<"\n";
			assert(j!=n);

			IntegerVectorList preCheck;
			for(IntegerVectorList::const_iterator l=F.begin();l!=F.end();l++)
				if((*k)[j]*(*l)[1]>0)
				{
		//			debug<<"YES\n";
					int LCM=(*k)[j]*(*l)[1]/gcdGFAN((*k)[j],(*l)[1]);
					int s1=LCM/(*k)[j];
					int s2=LCM/(*l)[1];
					if(s1<0)s1=-s1;
					if(s2<0)s2=-s2;
		//			debug<<"YES\n";
					IntegerVector v=s1*concatenation(concatenation(k->subvector(0,i),IntegerVector(1)),k->subvector(i,n-1))+s2*(*l)[0]*IntegerVector::standardVector(n,i);
					preCheck.push_back(v);
				}
			log2 debug<<"Combined:\n"<<preCheck<<"\n";

			IntegerVectorList preCheck2;
			for(IntegerVectorList::const_iterator i=preCheck.begin();i!=preCheck.end();i++)
			{
				PolynomialSet temp=initialForms(I,*i);
				if(!temp.containsMonomialGenerator()&&!initialForms(additionalPolys,*i).containsMonomialGenerator())
					preCheck2.push_back(*i);
			}

			log2 debug<<"Combined2:\n"<<preCheck2<<"\n";
			if(preCheck2.size()==1)
				ret.splice(ret.end(),preCheck2);
			else
				toCheck.splice(toCheck.end(),preCheck2);

			{
				IntegerVectorList::iterator kk=k;
				kk++;
				complicatedList.erase(k);
				k=kk;
			}
		}
			else k++;
	}
#if 0
	for(IntegerVectorList::const_iterator k=complicatedList.begin();k!=complicatedList.end();k++)
	{
//		debug<<I.getRing()<<"\n"<<I;
		int j;
		for(j=0;j<n-1;j++)
			if((*k)[j]!=0)break;
		debug<<"Checking:"<<*k<<" with chosen coordinate:"<<j<<"\n";
		assert(j!=n);
		int j2=j;
		if(j2>=i)j2++;
		list<int> chosenVariables;
		chosenVariables.push_back(i);
		chosenVariables.push_back(j2);
		IntegerVectorList A;
		A.push_back(IntegerVector::allOnes(n)-IntegerVector::standardVector(n,i)-IntegerVector::standardVector(n,j2));
		A.push_back(IntegerVector::allOnes(n));
		MatrixTermOrder T(A);
		PolynomialSet J=I;

		debug<<"Doing elimination\n";
		buchberger(&J,T,true);
		PolynomialRing newRing2(newRing.getField(),2);
		PolynomialSet P=J.polynomialRingIntersection(newRing2,&chosenVariables);
		debug<<"Done eliminating\n";
		assert(P.size()==1);

		debug<<"Doing polyhedral computation\n";
		//		debug<<"PFRONT:"<<P.front()<<P.front().exponents()<<"\n";
				//		PolyhedralFan FF=PolyhedralFan::normalFanOfNewtonPolytope(P.front());
		PolyhedralFan FF=PolyhedralFan::bergmanOfPrincipalIdeal(P.front());
		IntegerVectorList F;
		if(FF.dimensionOfLinealitySpace()==1)
			{
				assert(FF.conesBegin()!=FF.conesEnd());
				IntegerVectorList l=FF.conesBegin()->generatorsOfLinealitySpace();
				assert(l.size()==1);
				F.push_back(l.front());
				F.push_back(-l.front());
//				debug<<"LINgen"<<F<<"\n";
			}
		else
		{
//			debug<<"Calling getrays\n";
			F=FF.getRays();
		}

		debug<<"Potential lifts"<<F<<"\n";

		IntegerVectorList preCheck;
		for(IntegerVectorList::const_iterator l=F.begin();l!=F.end();l++)
			if((*k)[j]*(*l)[1]>0)
			{
	//			debug<<"YES\n";
				int LCM=(*k)[j]*(*l)[1]/gcdGFAN((*k)[j],(*l)[1]);
				int s1=LCM/(*k)[j];
				int s2=LCM/(*l)[1];
				if(s1<0)s1=-s1;
				if(s2<0)s2=-s2;
	//			debug<<"YES\n";
				IntegerVector v=s1*concatenation(concatenation(k->subvector(0,i),IntegerVector(1)),k->subvector(i,n-1))+s2*(*l)[0]*IntegerVector::standardVector(n,i);
				preCheck.push_back(v);
			}
		debug<<"Combined:\n"<<preCheck<<"\n";

		IntegerVectorList preCheck2;
		for(IntegerVectorList::const_iterator i=preCheck.begin();i!=preCheck.end();i++)
		{
			PolynomialSet temp=initialForms(I,*i);
			if(!temp.containsMonomialGenerator())
				preCheck2.push_back(*i);
		}

		debug<<"Combined2:\n"<<preCheck2<<"\n";
		if(preCheck.size()==1)//preCheck2!!!
			ret.splice(ret.end(),preCheck2);
		else
			toCheck.splice(toCheck.end(),preCheck2);
	}
#endif


	log2 debug<<"OldtoCheck:"<<toCheck;
	{
		IntegerVectorList toCheck2;
		for(IntegerVectorList::const_iterator i=toCheck.begin();i!=toCheck.end();i++)
		{
			PolynomialSet temp=initialForms(additionalPolys,*i);
			if(!temp.containsMonomialGenerator())
				toCheck2.push_back(*i);
		}
		toCheck=toCheck2;
	}


	log2 debug<<"Checking to check list\n";

	log2 debug<<"Already produced:"<<ret;
	log2 debug<<"toCheck:"<<toCheck;

	for(IntegerVectorList::const_iterator i=toCheck.begin();i!=toCheck.end();i++)
	{
		debug<<"Tropical variety membership test of ray:"<<*i<<"\n";
		stat_isInTropicalVariety++;
		if(isInTropicalVariety(I,*i,false))
		{
//			debug<<"PUSHING\n";
			ret.push_back(*i);
			if(earlyExit)return ret;
		}
		debug<<"Done membership test.\n";
	}



	/*	debug<<"FF"<<F<<"\n";
		debug<<"j2"<<j2<<"\n";
		debug<<"j"<<j<<"\n";
		debug<<"*k"<<*k<<"\n";
		//FF.printWithIndices(&debug);
*/
	log2
	{
		debug<<"At dimension "<<n<<"\n";
		debug<<"To Check\n"<<toCheck;
		debug<<"Complicated"<<complicatedList;

		debug<<"Summary n="<<n<<":\n";
		debug<<"isInTropicalVariety:"<<stat_isInTropicalVariety<<"\n";
		debug<<"#Rays for which getting finite lifts is easy:"<<stat_finiteLiftEasy<<"\n";
		debug<<"#Rays for which getting finite lifts is hard:"<<stat_finiteLiftHard<<"\n";


		debug<<"Returning"<<ret<<"\n";
	}

/*	{
		for(IntegerVectorList::const_iterator i=ret.begin();i!=ret.end();i++)
		{
			PolynomialSet iI=initialIdealNonHomogeneous(I,*i,true);
			PolynomialSet isat=nonHomogeneousSaturation(iI);
			debug<<*i<<" multiplicity: "<<multiplicity(isat)<<"\n";
		}
	}*/

	return ret;
}
