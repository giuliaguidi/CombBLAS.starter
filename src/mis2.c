#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

/* S = string matrix,
 * A = restriction matrix,
 * X = hop matrix,
 * nr1 = number of reads before contraction
 * nr2 = number of reads after contraction */

uint32_t *Sv, *Av, *Xv;
int      *Sr, *Sc, *Ar, *Ac, *Xr, *Xc;
int      nr1, nr2, Snnz, Annz, Xnnz;

void dense_to_csr(uint32_t **Mv, int **Mr, int **Mc, uint32_t *Mdense, int m, int n, int nnz);
void print_csr(uint32_t *Mv, int *Mr, int *Mc, int n);

int main(int argc, char *argv[])
{
    FILE     *fp;
    char     line[1024+1];
    int      i, j, k;
    uint32_t L, m, n, nnz, val;
    uint8_t  dir;

    assert(argc == 2);
    assert(fp = fopen(argv[1], "r"));

    do { fgets(line, 1024, fp); } while (*line == '%');

    sscanf(line, "%u %u %u\n", &m, &n, &nnz);
    assert(m == n);

    nr1 = (int)m;
    uint32_t *S = calloc(nr1*nr1, sizeof(*S));

    while (fscanf(fp, "%d %d %s %u\n", &i, &j, &dir, &L) == 4)
        S[--i*nr1 + --j] = (dir&3) | (L<<2);
        /* S[--i*nr1 + --j] = L; */

    fclose(fp);

    Snnz = (int)nnz;
    dense_to_csr(&Sv, &Sr, &Sc, S, nr1, nr1, nnz);
    free(S);

    print_csr(Sv, Sr, Sc, nr1);

    uint8_t *y = calloc(nr1, sizeof(*y));
    uint8_t *z = calloc(nr1, sizeof(*z));

    nr2 = 0;
    for (i = 0; i < nr1; ++i) {
        for (k = Sr[i]; k < Sr[i+1]; ++k) {
            if (Sv[k]) y[i]++;
        }
        if (y[i]==1) { nr2++; z[i] = 1; }
    }

    printf("\n\n");
    for (i = 0; i < nr1; ++i)
        printf("%d %d\n", i, z[i]);


    /* MIS2 here */

    return 0;
}

void print_csr(uint32_t *Mv, int *Mr, int *Mc, int n)
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int val = 0;
            for (int k = Mr[i]; k < Mr[i+1]; ++k) {
                if (Mc[k] == j) {
                    val = Mv[k];
                    break;
                }
            }
            printf("% 6d ", val);
        }
        printf("\n");
    }
}

void dense_to_csr(uint32_t **Mv, int **Mr, int **Mc, uint32_t *Mdense, int m, int n, int nnz)
{
    uint32_t *_Mv, val;
    int *_Mr, *_Mc, i, j, k;

    _Mv = calloc(nnz, sizeof(*_Mv));
    _Mc = calloc(nnz, sizeof(*_Mc));
    _Mr = calloc(m+1, sizeof(*_Mr));

    k = 0;
    for (i = 0; i < m; ++i) {
        _Mr[i] = k;
        for (j = 0; j < n; ++j) {
            if ((val = Mdense[i*n + j])) {
                _Mv[k] = val;
                _Mc[k] = j;
                ++k;
            }
        }
    }

    assert(k == nnz);
    _Mr[m] = k;

    *Mv = _Mv; *Mr = _Mr; *Mc = _Mc;
}
