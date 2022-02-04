import numpy as np
from numpy.linalg import svd
from scipy import linalg

#Helper function to obtain null space
def nullspace(A, atol=1e-5, rtol=0):
    A = np.atleast_2d(A)
    u, s, vh = svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    ss = s[nnz:]
    return ns, ss

dt = 0.000001
storage_time = 10000000
time = dt*storage_time
cells = 32

#build matrices
k_hit = np.zeros((cells,cells))
for i in range(cells):
    data = np.genfromtxt(str(i)+"/k_hits.txt")
    k_hit += data
k_hit /= time
k_hit = k_hit.T-np.diag(np.sum(k_hit,axis=1))
v,w = nullspace(k_hit)
index = np.argmin(w)
zl = np.real(v[:,index])
zl = zl/np.sum(zl)

width = 0.25
r_wca = 2**(1.0/6.0)
r_ext = r_wca+2*width
divide = int(cells/2)
prob_A = np.sum(zl[0:divide])
prob_B = np.sum(zl[divide:])
np.savetxt("zl.txt", zl)
np.savetxt("prob_A.txt", np.array((prob_A,)))
np.savetxt("prob_B.txt", np.array((prob_B,)))
