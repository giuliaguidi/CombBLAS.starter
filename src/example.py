import random
import numpy as np

S = np.zeros((8,8))
#  S[0,1] = 2
#  S[1,0] = 2
#  S[1,2] = 4
#  S[2,1] = 4
#  S[2,3] = 3
#  S[3,2] = 3
#  S[3,5] = 4
#  S[5,3] = 4
#  S[2,4] = 4
#  S[4,2] = 2
#  S[4,6] = 2
#  S[6,4] = 4
#  S[6,7] = 4
#  S[7,6] = 4

S[0,1] = 1
S[1,0] = 1
S[1,2] = 1
S[2,1] = 1
S[2,3] = 1
S[3,2] = 1
S[3,5] = 1
S[5,3] = 1
S[2,4] = 1
S[4,2] = 1
S[4,6] = 1
S[6,4] = 1
S[6,7] = 1
S[7,6] = 1

def MIS2(S):
    S2 = S@S + S
    n = S.shape[0]
    I = set([])
    C = set(range(n))
    while len(C) > 0:
        v = random.choice(list(C))
        N = set(np.where(S2[v,:]!=0)[0])
        N.add(v)
        C.difference_update(N)
        I.add(v)
    return list(I)

def contract(S):
    I = MIS2(S)
    m = S.shape[0]
    n = len(I)
    A = np.zeros((m,n))
    for i,v in enumerate(I):
        A[v,i] = 1
    A += S@A
    B = (A.T@S)@A
    return A.T, B




