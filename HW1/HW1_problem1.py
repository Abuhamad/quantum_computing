'''

Author: Mohammed Abuhamad
UCFID: 5013840
Description : 
			Assignment 1: Problem 1
			Compute the eigenbases of the Pauli operators 
			and check that they form mutually unbiased bases.

'''

import numpy as np
import numpy.linalg as la
np.set_printoptions(precision=2)

# pauli operators
sigma_x = np.array([[0, 1],[1, 0]])
sigma_y = np.array([[0, -1j],[1j, 0]])
sigma_z = np.array([[1, 0],[0, -1]])

#for these pauli matrices eigenvectors are eginbases
def get_bases(matrix):
    eigenvalues, eigenvectors = la.eig(matrix)
    return eigenvectors

#check whether two eigenbases are mutually unbiased (MUB)
def is_MUB (eigenbasis1, eigenbasis2):
    return np.around(np.square(la.norm(np.matmul(eigenbasis1, eigenbasis2))), decimals=1) == 1.0/len(eigenbasis1)

#get the eigenbases of matrices and check those bases
def check_bases_two_matrices (matrix1, matrix2):
    bases1 = get_bases(matrix1)
    bases2 = get_bases(matrix2)
    for basis1 in bases1:
        for basis2 in bases2:
            print( "Are eigenbases {} and {} mutually unbiased :> {}\n".format(basis1, basis2, is_MUB(basis1, basis2)))

if __name__ == '__main__':

	print ("Eigenbases of sigma_x are :{} and {} \n".format(get_bases(sigma_x)[:,0], get_bases(sigma_x)[:,1]))
	print ("Eigenbases of sigma_y are :{} and {} \n".format(get_bases(sigma_y)[:,0], get_bases(sigma_y)[:,1]))
	print ("Eigenbases of sigma_z are :{} and {} \n".format(get_bases(sigma_z)[:,0], get_bases(sigma_z)[:,1]))
	check_bases_two_matrices(sigma_x, sigma_y)
	check_bases_two_matrices(sigma_x, sigma_z)
	check_bases_two_matrices(sigma_y, sigma_z)