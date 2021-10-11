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

// TODO: create semiring

int main
{
    // TODO: we need to create B
    PSpMat<int>::MPI_DCCols BT = B; // build matrix
    
    /* how to create transpose matrix */
    BT.Transpose();
    
    /* how to apply a functor (semiring operation) to each nonzero in the matrix */
    BT.Apply(OverhangTSRing<int, int>()); 

    /* how to symmetricize a lower/upper triangular matrix */
    if(!(BT == B))
    {
        B += BT;
    }
  
    /* debug print */
    B.PrintInfo();
    
    /* print matrix to file in matrix market format */
    B.ParallelWriteMM("bridget-test.mm", true);     
	
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
                     
                   


