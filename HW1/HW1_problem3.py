# -*- coding: utf-8 -*-
'''

Author: Mohammed Abuhamad
UCFID: 5013840
Description : 
			Assignment 1: Problem 3 (Unitary error basis).

'''

import numpy as np
import numpy.linalg as la
np.set_printoptions(precision=4)


#define the matrices X and Z ∈ C^(d×d):
def define_matrices_X_Z(d):
	w = np.exp((2*180)/d) # primitive d-th root of unity
	X = np.zeros((d, d))
	Z = np.zeros((d, d))
	
	#construct X ∈ C^(d×d) = ∑|k+ 1〉〈k| from 0 to d-1
	for k in range (0, d):
		x1 = np.zeros(d)
		x2 = np.zeros(d)
		x1[(k+1)%d]= 1
		x2[k] = 1
		X = np.add(X, np.outer(x1.T, x2))

	#construct Z ∈ C^(d×d) = ∑ω^l|l〉〈l|
	for l in range (0, d):
		z1 = np.zeros(d)
		z2 = np.zeros(d)
		z1[l]= w**l
		z2[l] = 1
		Z = np.add(Z, np.outer(z1.T, z2))

	return X, Z

#construct M^(a,b) ∈ C^(d×d) = X^a Z^b for a,b ∈ {0,...,d−1}
def define_M (X, Z, a, b):
	return np.inner(la.matrix_power(X,a),la.matrix_power(Z,b))

#show that M^(a,b) from orthonormal basis with respect to the trace inner product
def check_basis(M):
	#get the bases:
	eigenvalues, eigenvectors = la.eig(M)
	
	'''
	#print all basis
	for idx, eigenvector in enumerate(eigenvectors):
		print('basis {}: {} \n'.format(idx, eigenvector))
	'''

	#check every pair (orthonormal if <v_i,v_j> = 0 when i!=j) and <v_i,v_i> = 1
	check = 0
	for idx in range(len(eigenvectors)-1):
		for jdx in range(idx, len(eigenvectors)):
			
			if idx == jdx:
				continue
				
				if np.around(np.inner(eigenvectors[idx].T, eigenvectors[jdx]), decimals = 1) == 1:
					check +=1
				else:
					print ('basis {} and {} are not orthonormal'.format(eigenvectors[idx], eigenvectors[jdx]))
				
			else:
				
				if np.around(np.inner(eigenvectors[idx].T, eigenvectors[jdx]), decimals = 1) == 0:
					check +=1

				else:
					print ('basis {} and {} are not orthonormal'.format(eigenvectors[idx], eigenvectors[jdx]))
	
	#if checks = len(eigenvectors) * (len(eigenvectors) - 1) / 2 then all pairs are checked
	if check == (len(eigenvectors) * (len(eigenvectors) - 1) / 2):
		print('All {} vector cobinations are orthonormal'.format(check))
	else:
		print ('Not all basis are orthonormal')




if __name__ == '__main__':
	#define the dimintion of basics matrices d
	# define a and b to construct M
	d = 3
	a = 1
	b = 1

	#get X and Z
	X, Z = define_matrices_X_Z(d)
	#get M
	M = define_M (X, Z, a, b)
	#check the bases of M
	check_basis(M)

