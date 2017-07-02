/*
 * gfanlib_tableau.h
 *
 *  Created on: May 5, 2017
 *      Author: anders
 */

#ifndef GFANLIB_TABLEAU_H_
#define GFANLIB_TABLEAU_H_

#include "gfanlib_matrix.h"
#include "gfanlib_circuittableint.h"

/*
 * Simplex algorithm without objective function.
 * The original matrix does not have to be stored.
 * At any time a B denotes a submatrix of columns
 * indexed by basisIndices. This matrix also does
 * not have to be stored. The adjoint matrix adj(B)
 * is also unknown, but its product with the original
 * matrix is kept updated as the columns of the basis
 * B are exchanged. Schrijver version of Bland's rule is used.
 */

namespace gfan{

template <class matrixType> std::string matrixToString(matrixType const &m)
{
	int h=m.getHeight();
	int w=m.getWidth();
	vector<vector<string> > theStrings;
	for(int i=0;i<h;i++)
	{
		vector<string> row;
		for(int j=0;j<w;j++)row.push_back(m[i][j].toString());
		theStrings.push_back(row);
	}
	vector<int> columnWidths(w);
	for(int j=0;j<w;j++)
	{
		int width=1;
		for(int i=0;i<h;i++)if(width<theStrings[i][j].size()+1)width=theStrings[i][j].size()+1;
		columnWidths[j]=width;
	}
	std::stringstream s;
	for(int i=0;i<h;i++)
	{
		for(int j=0;j<w;j++)
		{
			for(int k=theStrings[i][j].size();k<columnWidths[j];k++)s<<" ";
			s<<theStrings[i][j];
		}
		s<<"\n";
	}
	return s.str();
}

template <class mvtyp> class Tableau{
	public:
		Matrix<mvtyp> combinedMatrix; // Adj(B)*M
		std::vector<int> basisIndices;
		std::vector<bool> inBasis;
		mvtyp determinantOfBasis;
		Matrix<mvtyp> original; // Delete
		Tableau(Matrix<mvtyp> const &M, bool appendIdentity, bool appendAllOnes=false):
			basisIndices(M.getHeight()+appendAllOnes)
		{
			if(appendAllOnes)
			{
				Matrix<mvtyp> M2=M;
				M2.appendRow(Vector<mvtyp>::allOnes(M2.getWidth()));
				combinedMatrix=combineLeftRight(M2.identity(M.getHeight()+1),M2);
			}
			else
				combinedMatrix=combineLeftRight(M.identity(M.getHeight()),M);
			inBasis=std::vector<bool>(combinedMatrix.getWidth());
			for(int i=0;i<M.getHeight()+appendAllOnes;i++){basisIndices[i]=i;inBasis[i]=true;}
			determinantOfBasis=1;
			original=combinedMatrix;
		}
		static std::string vectorToString(std::vector<int> v)
		{
			std::stringstream s;
			s<<"(";
			for(int i=0;i<v.size();i++)
			{
				if(i)s<<",";
				s<<v[i];
			}
			s<<")";
			return s.str();
		}
		std::string toString()
		{
			return "Tableau:\n"+
//			combinedMatrix.toString()+
			matrixToString<Matrix<mvtyp> >(combinedMatrix)+
			"Basis:\n"+vectorToString(basisIndices)+"\nDeterminant:\n"+determinantOfBasis.toString()+"\n";
		}
		void check()
		{
			int d=original.getHeight();
			for(int i=0;i<d;i++)
				for(int j=0;j<d;j++)
				{
					mvtyp s=0;
					for(int k=0;k<d;k++)
						s+=original[k][basisIndices[i]]*combinedMatrix[j][k];
/*					if(i==j)
						std::cerr<<"Equ:"<<s.toString()<<"\n";
					else
						std::cerr<<"zero:"<<s.toString()<<"\n";*/
				}

		}
		void exchange(int i, int j) // make assigment basis[i]=j, updating matrix
		{
			std::cerr<<"EXCHANGE("<<i<<","<<j<<")\n"<<toString();
			check();
			mvtyp detMultiplier=combinedMatrix[i][j];
			std::cerr<<"MUL"<<detMultiplier<<"\n";
			for(int k=0;k<combinedMatrix.getHeight();k++)
//				combinedMatrix.madd(i,-combinedMatrix[k][j],k);
				if(k!=i)
			{
				mvtyp temp=-combinedMatrix[k][j];
				debug<<"Madd("<<i<<","<<temp.toString()<<","<<k<<"\n";
				for(int a=0;a<combinedMatrix.getWidth();a++)
					{
//					combinedMatrix[k][a]+=temp*combinedMatrix[k][a]/determinantOfBasis;
//					combinedMatrix[k][a]+=dotDiv(0,0,temp,combinedMatrix[i][a],determinantOfBasis);

//					std::cerr<<combinedMatrix[i][j]<<combinedMatrix[k][a]<<temp<<combinedMatrix[i][a]<<determinantOfBasis<<"\n";
//					std::cerr<<int32_t((((int32_t)(((uint64_t)(((int64_t)combinedMatrix[i][j].v)*((int64_t)combinedMatrix[k][a].v)+((int64_t)temp.v)*((int64_t)combinedMatrix[i][a].v)))>>CircuitTableInt32::Divisor(determinantOfBasis).shift))))<< "sh"<<CircuitTableInt32::Divisor(determinantOfBasis).shift<<"inv"<<	CircuitTableInt32::Divisor(determinantOfBasis).multiplicativeInverse<<"\n";
//					std::cerr<<int32_t((((int32_t)(((uint64_t)(((int64_t)combinedMatrix[i][j].v)*((int64_t)combinedMatrix[k][a].v)+((int64_t)temp.v)*((int64_t)combinedMatrix[i][a].v)))>>CircuitTableInt32::Divisor(determinantOfBasis).shift)))*CircuitTableInt32::Divisor(determinantOfBasis).multiplicativeInverse)<<"\n";
							combinedMatrix[k][a]=dotDiv(combinedMatrix[i][j],combinedMatrix[k][a],temp,combinedMatrix[i][a],determinantOfBasis/*1??*/);
//					std::cerr<<combinedMatrix[k][a]<<"\n";
					}
			}
			{
				mvtyp temp=1;//combinedMatrix[i][j];
				for(int a=0;a<combinedMatrix.getWidth();a++)
					combinedMatrix[i][a]=temp*combinedMatrix[i][a];
			}
			inBasis[basisIndices[i]]=false;
			basisIndices[i]=j;
			inBasis[j]=true;
			determinantOfBasis=detMultiplier;
		}
	};
	template<class mvtyp> class TableauSolver:public Tableau<CircuitTableInt32>{
	public:
		vector<bool> ignoreRows;
		vector<bool> nonExtreme;
		vector<bool> inLineality;
//		vector<bool> permanentBasisMember;
		bool appendIdentity;
		std::string toString()
		{
			std::stringstream s;
			s<<"nonExtreme:";for(int i=0;i<combinedMatrix.getWidth();i++)s<<nonExtreme[i];s<<"\n";
			s<<"inLinealin:";for(int i=0;i<combinedMatrix.getWidth();i++)s<<inLineality[i];s<<"\n";
//			s<<"permanentB:";for(int i=0;i<combinedMatrix.getWidth();i++)s<<permanentBasisMember[i];s<<"\n";
			s<<"ignoreRows:";for(int i=0;i<combinedMatrix.getHeight();i++)s<<ignoreRows[i];s<<"\n";
			return Tableau<CircuitTableInt32>::toString()+s.str();
//			return "Tableau:\n"+combinedMatrix.toString()+"Basis:\n"+vectorToString(basisIndices)+"\nDeterminant:\n"+determinantOfBasis.toString()+"\n";
		}
		TableauSolver(Matrix<mvtyp> const &M, bool appendIdentity_):
			Tableau<CircuitTableInt32>(M, appendIdentity_,true),
			ignoreRows(M.getHeight()),
			nonExtreme(M.getHeight()+M.getWidth()+1),
			appendIdentity(appendIdentity_),
			inLineality(M.getHeight()+M.getWidth()+1)//,
	//		permanentBasisMember(M.getHeight()+M.getWidth()+1)
		{
		}
		bool loose(int basisIndex, int lowerIndex, int upperIndex, vector<bool> &nonExtreme)
		{
			int basisIndex2=-1;
			for(int i=0;i<basisIndices.size();i++)if(basisIndex==basisIndices[i]){basisIndex2=i;}

			std::cerr<<"loose called"<<basisIndex<<"lower"<<lowerIndex<<"upper"<<upperIndex<<"s\n";
			for(int i=lowerIndex;i<upperIndex;i++)
				if(!nonExtreme[i]&&!inLineality[i])
				if((!inBasis[i])&&combinedMatrix[basisIndex2][i].isNonZero())
				{
					exchange(basisIndex2,i);
					return true;
				}
			return false;
		}
		bool towards(int i, int lowerIndex, int upperIndex, int ignore, vector<bool> const &ignoreRows, vector<bool> &nonExtreme,bool ignoreAddedRow)
		{
			//i is not supposed to be in current basis
			//return true when simplicial cone containing i is found
			assert(inBasis[i]==false);
			int targetIndex=i;
			while(1)
			{
				//check whether targetIndex column is in simplex of basis
				int violatedIndex2=1000000000;
				int violatedIndex=-1;
				std::cerr<<targetIndex<<"\n";
				for(int i=0;i<combinedMatrix.getHeight()-ignoreAddedRow;i++)
					if(!ignoreRows[i])
				{
					if(determinantOfBasis.isNegative())
						if(combinedMatrix[i][targetIndex].isPositive()){if(basisIndices[i]<violatedIndex2){violatedIndex=i;violatedIndex2=basisIndices[i];}}
					if(determinantOfBasis.isPositive())
						if(combinedMatrix[i][targetIndex].isNegative()){if(basisIndices[i]<violatedIndex2){violatedIndex=i;violatedIndex2=basisIndices[i];}}
				}
				std::cerr<<violatedIndex<<"\n";
				if(violatedIndex==-1)
				{
					std::cerr<<"Inside cone!\n";
//					nonExtreme[i]=true;
					return true;
				}
				else
				{
					std::cerr<<"Violated:"<<violatedIndex2<<"\n";
					int i;
					for(i=lowerIndex;i<upperIndex;i++)
						if(!nonExtreme[i])
							if(!inLineality[i])
						if(i!=targetIndex)
						if(!inBasis[i])
					{
						if(combinedMatrix[violatedIndex][i].isNegative()&&determinantOfBasis.isPositive() ||
								combinedMatrix[violatedIndex][i].isPositive()&&determinantOfBasis.isNegative()
								)
						{
							exchange(violatedIndex,i);
							std::cerr<<toString();
							break;
						}
					}
					if(i==upperIndex)
					{
						std::cerr<<"Separating hyperplane found:\n";
						std::cerr<<combinedMatrix.submatrix(violatedIndex,0,violatedIndex+1,combinedMatrix.getHeight()).toString()<<"\n";
						return false;
					}
				}
			}
			return false;
		}
		void coneInfo()//test routine for finding information about a cone given by inequalities.
		{
			std::cerr<<toString();

			int lowerIndex=combinedMatrix.getHeight();
			int upperIndex=combinedMatrix.getWidth();

			std::cerr<<"WE FIRST FIGURE OUT THE ORTHORGONAL COMPLEMENT.\n";
			int d=combinedMatrix.getHeight();
			for(int i=0;i<d-1;i++)
			{
				int j;
				for(j=lowerIndex;j<upperIndex;j++)
					if(combinedMatrix[i][j].isNonZero())
						break;
				if(j==upperIndex)
				{
					std::cerr<<"FLAT MUST IGNORE ROW"<<i<<"\n";
					std::cerr<<"CAUSE:"<<combinedMatrix.submatrix(i,0,i+1,d).toString()<<"\n";
					ignoreRows[i]=true;
				}
				else
				{
					this->exchange(i,j);
				}
				std::cerr<<toString();
				std::cerr<<"inBasis";
				for(int j=0;j<inBasis.size();j++)std::cerr<<inBasis[j];
			}
			std::cerr<<"NOW WE MUST FIND LINEALITY SPACE.\n";
			do
			{
				//we need to check if some is in span of added vector
				for(int i=lowerIndex;i<upperIndex;i++)
				{
					bool isZeroColumn=true;
					for(int j=0;j<combinedMatrix.getHeight()-1;j++)
						if(!ignoreRows[j])if(combinedMatrix[j][i].isNonZero())isZeroColumn=false;
					if(isZeroColumn)
					{
						std::cerr<<"InLinelity"<<i<<"\n";
						inLineality[i]=true;
					}
				}
				if(loose(d-1,lowerIndex,upperIndex,nonExtreme))
				{
					std::cerr<<"WE LOST ADDED\n";
					bool foundSimplex=towards(d-1,lowerIndex,upperIndex,d-1,ignoreRows,nonExtreme,false);
					std::cerr<<"FS"<<foundSimplex<<"\n";

					if(!foundSimplex)break;
					std::cerr<<toString();
					// Now we have found a positive circuit involving the basis elements with non-zero coefficient in special column and the special column
					// We must mark all that are in span
					for(int i=lowerIndex;i<upperIndex;i++)
					{
						bool isInSpan=true;
						for(int j=0;j<d;j++)
							if(!ignoreRows[j])
								if(combinedMatrix[j][d-1].isZero())
									if(combinedMatrix[j][i].isNonZero())isInSpan=false;
						if(isInSpan)inLineality[i]=true;
					}
					// We must swap in the special column
					//check that we can just use normal position. If not things will probably need to get complicated
					assert(combinedMatrix[d-1][d-1].isNonZero());
					int makePermanentIndex;
					for(makePermanentIndex=0;makePermanentIndex<d;makePermanentIndex)if(combinedMatrix[makePermanentIndex][d-1].isNonZero())break;
					std::cerr<<toString();
					std::cerr<<"EXCHANGE!\n";
					exchange(d-1,d-1);
					std::cerr<<toString();

					// We must make one of the other circuit members a permanent member of the basis by setting ignoreRow
//					assert(0);
					assert(makePermanentIndex!=d-1);
					assert(makePermanentIndex!=d);
					ignoreRows[makePermanentIndex]=true;
				}
				else
				{
					std::cerr<<"WHAT DO WE DO NOW?\n";
					assert(0);
				}
			}while(1);


			assert(combinedMatrix[d-1][d-1].isNonZero());
			std::cerr<<toString();
			std::cerr<<"MAKE SURE THAT d-1 is in basis!\n";
			exchange(d-1,d-1);
			std::cerr<<toString();

			std::cerr<<"NOW WE FIND EXTREME RAYS.\n";


			for(int i=0;i<d;i++)nonExtreme[i]=true;
			for(int i=lowerIndex;i<upperIndex;i++)if(inLineality[i])nonExtreme[i]=true;


			for(int i=lowerIndex;i<upperIndex;i++)
				if(!inLineality[i])
			{
//				std::cerr<<"A1\n";
				//we wish to check if ith column is extreme
				std::cerr<<	"index:"<<i<<"\n";
				std::cerr<<"inBasis";
				for(int j=0;j<inBasis.size();j++)std::cerr<<inBasis[j];
				std::cerr<<this->toString();
				if(inBasis[i] && !loose(i,lowerIndex,upperIndex,nonExtreme))//we first must remove it from the basis to check
				{
					std::cerr<<this->toString();
//					std::cerr<<"A2\n";
					std::cerr<<"Is extreme because it cannot be replaced\n";
					nonExtreme[i]=false;
					//check scaling?
				}
				else
				{
					std::cerr<<this->toString();
//					std::cerr<<"A3\n";
//					debug<<int(
					bool foundSimplex=towards(i,lowerIndex,upperIndex,i,ignoreRows,nonExtreme,true);//<<"\n";
					if(foundSimplex)nonExtreme[i]=true;
				}

			}
//			std::cerr<<"A4\n";

			// Is the cone pointed?
/*			for(int i=0;i<combinedMatrix.getHeight();i++)
			{
				if(combinedMatrix[i][targetIndex].isPositive()){if(basisIndices[i]<violatedIndex2){violatedIndex=i;violatedIndex2=basisIndices[i];}}
			}
	*/
			std::cerr<<"Extreme:\n";
			for(int i=0;i<nonExtreme.size();i++)
				if(!nonExtreme[i])std::cerr<<i<<"\n";
		}
		void solve()
		{
			// try to express
			int lowerIndex=combinedMatrix.getHeight();
			int upperIndex=combinedMatrix.getWidth()-1;
			int targetIndex=upperIndex;
			while(1)
			{
				//check whether last column is in simplex of basis
				int violatedIndex2=1000000000;
				int violatedIndex=-1;
				std::cerr<<targetIndex<<"\n";
				for(int i=0;i<combinedMatrix.getHeight();i++)
				{
					if(determinantOfBasis.isNegative())
						if(combinedMatrix[i][targetIndex].isPositive()){if(basisIndices[i]<violatedIndex2){violatedIndex=i;violatedIndex2=basisIndices[i];}}
					if(determinantOfBasis.isPositive())
						if(combinedMatrix[i][targetIndex].isNegative()){if(basisIndices[i]<violatedIndex2){violatedIndex=i;violatedIndex2=basisIndices[i];}}
				}
				std::cerr<<violatedIndex<<"\n";
				if(violatedIndex==-1)
				{
					std::cerr<<"Inside cone!\n";
					return;
				}
				else
				{
					std::cerr<<"Violated:"<<violatedIndex2<<"\n";
					int i;
					for(i=lowerIndex;i<upperIndex;i++)
						if(!inBasis[i])
					{
						if(combinedMatrix[violatedIndex][i].isNegative()&&determinantOfBasis.isPositive() ||
								combinedMatrix[violatedIndex][i].isPositive()&&determinantOfBasis.isNegative()
								)
						{
							exchange(violatedIndex,i);
							std::cerr<<toString();
							break;
						}
					}
					if(i==upperIndex)
					{
						std::cerr<<"Inside cone!\n:Separating hyperplane found:\n";
						return;
					}
				}
			}
		}
	};
}

#endif /* GFANLIB_TABLEAU_H_ */
