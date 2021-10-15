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

#include <assert.h>
#include <stdio.h>

// template <class NT>
// class PSpMat
// {
// public:
  // typedef combblas::SpTuples <int64_t, NT> Tuples;
  // typedef combblas::SpDCCols <int64_t, NT> DCCols;
  // typedef combblas::SpParMat <int64_t, NT, DCCols> MPI_DCCols;
  // typedef std::tuple<int64_t, int64_t, NT *> ref_tuples;
// };

int main(int argc, char **argv)
{

    /* S = string matrix after transitive reduction */

    /* ST = S.Transpose() */

    /* ST.Appy(ChangeEntrySR<>) //  check what this means */

    /* do
     *     // CCGrid& CMG is related to the 3d SpGEMM
     *     // S is the string matrix, A and AT represent the "Assembly Matrix"
     *     // RestrictionOp take sthe S matrix as input, finds MIS2, and
     *     // contracts it, deep copying it into A/AT.
     *     RestrictionOp(S, A, AT)
     * while (only edges that remain are edges at some branching point)
     *
     * // Once we stop contracting, want to align sequences of long aggregates
     * // starting from branching point.
     */

    int nproc, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    FILE *fp;

    assert(argc == 2);
    assert(fp = fopen(argv[1], "r"));

    char line[1024 + 1];
    int64_t i, j;
    uint32_t L, m, n, nnz;
    uint8_t dir;

    typedef std::pair<uint8_t, uint32_t> > ReadEdge;

    std::vector<int64_t> cols, rows;
    std::vector<ReadEdge> vals;

    sscanf(line, "%u %u %u\n", &m, &n, &nnz);
    assert(m == n);

    while (fscanf(fp, "%d %d %u %u\n", &i, &j, &dir, &L) == 4) {
        rows.push_back(i);
        cols.push_back(j);
        vals.push_back(std::make_pair(dir, L));
    }

    fclose(fp);

    std::shared_ptr<combblas::CommGrid> fullWorld;
    fullWorld.reset(new combblas::CommGrid(MPI_COMM_WORLD, 0, 0));

    combblas::FullyDistVec<int64_t, int64_t> drows(rows, fullWorld);
    combblas::FullyDistVec<int64_t, int64_t> dcols(cols, fullWorld);
    combblas::FullyDistVec<int64_t, std::pair<uint8_t, uint32_t> > dvals(vals, fullWorld);

    typedef SpParMat<int64_t, ReadEdge, combblas::SpDCCols<int64_t, ReadEdge> > StrMat;



    MPI_Finalize();
    return 0;
}
