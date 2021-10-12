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

#define MYCONSTANT 3 // @Bridget, this is just an example
#define MAXK 10      // @Bridget, change this based on your input matrix

//////////////////////////////////////////////////////////////////////////////////////
// DEFINITIONS AND SEMIRINGS                                                        // 
//////////////////////////////////////////////////////////////////////////////////////

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

const uint infplus(const int& a, const int& b) 
{
	uint inf = std::numeric_limits<uint>::max();
    if (a == inf || b == inf)
	{
    	return inf;
    }
    return a + b;
}

/* MinPlus semiring: given two parallel routes in the matrix, keep the shortest */
template <class T1, class T2, class OUT>
struct BridgetMinPlusSRing
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
	// static OUT axpy(T1 a, const T2& x, OUT& y)
	// {   
	// 	y = add(y, multiply(a, x));
	// }
};

template <class T>
struct BridgetZeroSR : unary_function <T, bool>
{
    bool operator() (const T& x) const { if(x == 0) return true; else return false; }
};

template <class T1, class T2, class OUT>
struct BridgetReduceMaxSRing : binary_function <T1, T2, OUT>
{
    OUT operator() (const T1& x, const T2& y) const
    {
        if(y > x) return static_cast<OUT>(y);
        else return static_cast<OUT>(x);
    }
};

template <class T, class OUT>
struct BridgetPlusConstSRing : unary_function <T, OUT>
{
    OUT operator() (T& x) const
    {
        return static_cast<OUT>(x + MYCONSTANT);
    }
};

/*! type definitions */
typedef BridgetMinPlusSRing <int, int, int> MinPlusSR_t;
typedef BridgetReduceMaxSRing <int, int, int> ReduceMSR_t;

//////////////////////////////////////////////////////////////////////////////////////
// MAIN FUNCITON                                                                    // 
//////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
	/*! Init MPI */
	int nproc;
	int rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	std::vector<int64_t> lcols, lrows;  // @Bridget, int64_t is the index type of your matrices right now (see PSpMat class definition above)
	std::vector<int> lvals;				// @Bridget, working on ints in this example

	/*! Create fictious matrix */
	for(int i = 0; i < MAXK; ++i)
    {
        for(int j = i; j < MAXK; ++j) // @Bridget, I'm filling the lower triangul matrix with 1s
        {
            lcols.push_back(j); // column ids
            lrows.push_back(i); // rows ids
            lvals.push_back(1); // values
        }
    }

	/*! Define the proc grid */
	shared_ptr<CommGrid> fullWorld;
	fullWorld.reset(new CommGrid(MPI_COMM_WORLD, 0, 0));

	/*! Create distributed sparse matrix */
    FullyDistVec<int64_t, int64_t> drows(lrows, fullWorld);
    FullyDistVec<int64_t, int64_t> dcols(lcols, fullWorld);
    FullyDistVec<int64_t, int> dvals(lvals, fullWorld);

    PSpMat<int>::MPI_DCCols A(MAXK, MAXK, drows, dcols, dvals, false);

#ifdef SYMM
	/*! Create a copy */
	PSpMat<int>::MPI_DCCols BT = B; // build matrix

	/*! How to create transpose matrix */
	BT.Transpose();

	/*! How to apply a functor (semiring operation) to each nonzero in the matrix */
	BT.Apply(OverhangTSRing<int, int>()); 

	/*! How to symmetricize a lower/upper triangular matrix */
	if(!(BT == B))
	{
	 += BT;
	}
#endif SYMM

	/*! Debug print */
	A.PrintInfo();
	
	/*! Print matrix to file in matrix market format */
	A.ParallelWriteMM("bridget-test-input.mm", true);     
	
	/*! Get some timing */
	uint nnz, prev;
	double timeA2 = 0, timeC = 0, timeI = 0, timeA = 0;

	/*! SpGEMM */
	double start = MPI_Wtime();
	
	PSpMat<int>::MPI_DCCols A2 = A;
	PSpMat<int>::MPI_DCCols B = Mult_AnXBn_DoubleBuff<MinPlusSR_t, int, PSpMat<int>::DCCols>(A, A2);
	
	timeA2 += MPI_Wtime() - start;

	/*! Prune evaluates a function, if true, remove that nonzero */ 
	B.Prune(BridgetZeroSR<int>(), true);

	start = MPI_Wtime();
	
	/*! It's creating an empty distributed vector mapped to the same processor grid as the matrix B */
	FullyDistVec<int64_t, int> Vec1(A.getcommgrid());

	int id = 0; 
	/*! This applies a reduction operation over rows of B and it returns a vector containing the result of the function evaluation */
	Vec1 = A.Reduce(Row, ReduceMSR_t(), id);

	/*! How to apply a functor (semiring operation) to each nonzero in the vector */
	Vec1.Apply(BridgetPlusConstSRing<int, int>());

	// @Bridget, it compiles/works fine until this point —I'm gonna leave the rest to you, feel free to reach out whenever needed.

	// /*! Basically the opposite od Reduce over Row in a matrix */
	// B1.DimApply(Row, Vec1, Bind2ndSR_t());
	
	// timeC += MPI_Wtime() - start;

	// start = MPI_Wtime();
	// bool isLogicalNot = false;
	
	// /*! Ewiseapply takes two matrixes as input and apply a function to their entries element wise */
	// PSpMat<bool>::MPI_DCCols I = EWiseApply<bool, PSpMat<bool>::DCCols>(B1, C, GreaterBinaryOp<int, int>(), isLogicalNot, id);

	// I.Prune(ZeroUnaryOp<bool>(), true); 
	
	// timeI += MPI_Wtime() - start;
	
	// start = MPI_Wtime();
	// isLogicalNot = true;
	
	// /*! Ewiseapply takes two matrixes as input and apply a function to their entries element wise */
	// B = EWiseApply<int, PSpMat<int>::DCCols>(B, I, EWiseMulOp<int, bool>(), isLogicalNot, true);

	// /*! Prune zero-valued overhang */
	// B.Prune(ZeroOverhangSR<int>(), true);
	// timeA += MPI_Wtime() - start;

	// B.PrintInfo();
	// B.ParallelWriteMM("bridget-test-output.mm", true); 

 #ifdef TIME
	
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

	/*! Finalize MPI */
	MPI_Finalize();

}
                     

