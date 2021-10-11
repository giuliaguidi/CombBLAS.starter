#include "CombBLAS/CombBLAS.h"

#include <sys/time.h> 
#include <iostream>
#include <functional>
#include <algorithm>
#include <vector>
#include <cmath>

#include <map>
#include <fstream>

#include <string>
#include <sstream>

/*! Namespace declarations */
using namespace combblas;
using namespace std;

// #define SYMM
// #define TIME

template <class T, std::size_t N>
ostream& operator<<(ostream& o, const array<T, N>& arr)
{
    copy(arr.cbegin(), arr.cend(), ostream_iterator<T>(o, " "));
    return o;
}

template <class NT>
class PSpMat
{
public:
  typedef combblas::SpTuples <int64_t, NT> Tuples;	
  typedef combblas::SpDCCols <int64_t, NT> DCCols;
  typedef combblas::SpParMat <int64_t, NT, DCCols> MPI_DCCols;
  typedef std::tuple<int64_t, int64_t, NT *> ref_tuples;
};

/* MinPlus semiring: given two parallel routes in the matrix, keep the shortest */
template <class T1, class T2, class OUT>
struct MinPlusSRing
{
	static OUT id() 	        { return std::numeric_limits<OUT>::max(); };
	static bool returnedSAID() 	{ return false; 	}
	static MPI_Op mpi_op() 		{ return MPI_MIN; 	};

	/* instead of add, return the min */
	static OUT add(const OUT& arg1, const OUT& arg2)
	{
		return min(arg1, arg2);
	}
	/* instead of multiply, return the sum */
	static OUT multiply(const T1& arg1, const T2& arg2)
	{
        	OUT sum = infplus(arg1, arg2);
           	return sum;
        } 
        else return id();
	}
	static void axpy(T1 a, const T2 & x, OUT & y)
	{   
		y = add(y, multiply(a, x));
	}
};

/*! type definitions */
typedef MinPlusBiSRing <int, int, int> MinPlusSR_t;

int main
{
	// TODO: we need to create B

#ifdef SYMM
	/* create a copy */
	PSpMat<int>::MPI_DCCols BT = B; // build matrix

	/* how to create transpose matrix */
	BT.Transpose();

	/* how to apply a functor (semiring operation) to each nonzero in the matrix */
	BT.Apply(OverhangTSRing<int, int>()); 

	/* how to symmetricize a lower/upper triangular matrix */
	if(!(BT == B))
	{
	 += BT;
	}
#endif SYMM

	/* debug print */
	B.PrintInfo();
	
	/* print matrix to file in matrix market format */
	B.ParallelWriteMM("bridget-test-input.mm", true);     
	
	/* variables to get some timing */
	uint nnz, prev;
	double timeA2 = 0, timeC = 0, timeI = 0, timeA = 0;

	/* mat mul
	 * C = B^2
	 */
	double start = MPI_Wtime();
	
	PSpMat<int>::MPI_DCCols B2 = B;
	PSpMat<int>::MPI_DCCols C = Mult_AnXBn_DoubleBuff<MinPlusSR_t, int, PSpMat<int>::DCCols>(B, B2);
	
	timeA2 += MPI_Wtime() - start;

	/* Prune evaluates a function, if true, remove that nonzero */ 
	C.Prune(ZeroOverhangSR<int>(), true);

	start = MPI_Wtime();
	
	/* It's creating an empty distributed vector mapped to the same processor grid as the matrix B */
	FullyDistVec<int64_t, int> Vec1(B.getcommgrid());

	int id = 0; 
	
	/* This applies a reduction operation over rows of B and it returns a vector containing the result of the function evaluation */
	Vec1 = B.Reduce(Row, ReduceMSR_t(), id);
	
	/* how to apply a functor (semiring operation) to each nonzero in the vector */
	Vec1.Apply(PlusFBiSRing<int, int>());

	/* basically the opposite od Reduce over Row in a matrix */
	B1.DimApply(Row, Vec1, Bind2ndSR_t());
	
	timeC += MPI_Wtime() - start;

	start = MPI_Wtime();
	bool isLogicalNot = false;
	
	/* ewiseapply takes two matrixes as input and apply a function to their entries element wise */
	PSpMat<bool>::MPI_DCCols I = EWiseApply<bool, PSpMat<bool>::DCCols>(B1, C, GreaterBinaryOp<int, int>(), isLogicalNot, id);

	I.Prune(ZeroUnaryOp<bool>(), true); 
	
	timeI += MPI_Wtime() - start;
	
	start = MPI_Wtime();
	isLogicalNot = true;
	
	/* ewiseapply takes two matrixes as input and apply a function to their entries element wise */
	B = EWiseApply<int, PSpMat<int>::DCCols>(B, I, EWiseMulOp<int, bool>(), isLogicalNot, true);

	/* Prune zero-valued overhang */
	B.Prune(ZeroOverhangSR<int>(), true);
	timeA += MPI_Wtime() - start;

	B.PrintInfo();

 #ifdef TIME
	B.ParallelWriteMM("bridget-test-output.mm", true); 
	
	double maxtimeA2, maxtimeC, maxtimeI, maxtimeA;
	
	MPI_Reduce(&timeA2, &maxtimeA2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&timeC,  &maxtimeC,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&timeI,  &maxtimeI,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&timeA,  &maxtimeA,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    	if(myrank == 0)
    	{
     		std::cout << "TransitiveReduction:TimeA2 = " << maxtimeA2 << std::endl;
      		std::cout << "TransitiveReduction:TimeC  = " <<  maxtimeC << std::endl;
      		std::cout << "TransitiveReduction:TimeI  = " <<  maxtimeI << std::endl;
      		std::cout << "TransitiveReduction:TimeA  = " <<  maxtimeA << std::endl;
    	}
#endif
}
                     
