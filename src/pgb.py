from random import random
from pygraphblas import *
from pygraphblas.types import Type, binop
from scipy.sparse import csr_matrix
import numpy as np

#  @select_op(FP64, FP64)
#  def random_le(i, x, v):
    #  if random() < v:
        #  return True
    #  return False

#  A = Matrix.sparse(FP64, 8, 8)
#  A[0,1] = 1
#  A[1,2] = 1
#  A[2,3] = 1
#  A[3,5] = 1
#  A[2,4] = 1
#  A[4,6] = 1
#  A[6,7] = 1
#  A += A.T

#  degree = A.reduce_vector(accum=FP64.PLUS)
#  prob = 1./(2*degree)
#  select = Vector.dense(BOOL, 8)

#  for i, v in enumerate(np.random.uniform(0, 1, 8)):
    #  select[i] = 1 if v < prob[i] else 0

#  neighbors = select & (A@select)

A = np.zeros((8, 8))
A[0,1] = 1
A[1,2] = 1
A[2,3] = 1
A[3,5] = 1
A[2,4] = 1
A[4,6] = 1
A[6,7] = 1
A += A.T

prob = 1/(2*A.sum(axis=1))
select = 1 * (np.random.uniform(0,1,8) < prob)

#  class StringEdge(Type):
    #  _base_name = "UDT"
    #  _numpy_t = None
    #  members = ["uint32_t L", "uint8_t dir"]

#  def read_strgraph(filename):
    #  with open(filename, "r") as f:
        #  next(f)
        #  n,_,nnz = (int(v) for v in next(f).rstrip().split())
        #  A = np.zeros((n,n), dtype=np.uint32)
        #  for k in range(nnz):
            #  i,j,d,L = (int(v) for v in next(f).rstrip().split())
            #  A[i-1,j-1] = (d&3)|(L<<2)
    #  return csr_matrix(A)

