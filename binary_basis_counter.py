#this script computes N_n, the number of distinct bases for R^n using only vectors whose entries are 0's and 1's for the n=2,3,4,5 cases
#to find N(n), we add the order of GL_n(2) to the epsilon term where the epsilon term is the number of nxn matrices of 0s and 1s with even determinant
#that is, N(n) = GL_n(2) + epsilon_term(n)
#there exists a closed-form solution for GL_n(q) for prime q. thus if we can find a closed-form solution for epsilon_term(n), we will have a closed-form solution for N(n)

import numpy as np
import math
import random
import itertools


def same_array_check(a,b):
    #input arrays a and b of same length
    #outputs true if the arrays are identical and false if not
    count = 0
    for i in range(0,len(a)):
        if a[i] == b[i]:
            count = count + 1
    if count == len(a):
        return True
    else:
        return False



n = 4 #dimension of R^n
data = np.ones((n,pow(2,n)-1)) #this is a matrix whose columns will represent the possible basis elements
count_epsilon = 0 #this will end up being the value of the epsion
numb_lin_indep_mat = 0 #this will end up being N(n)


#compute the order of GL_n(F_2)
order_GL = 1
for i in range(0,n):
    order_GL = order_GL * (pow(2,n)-pow(2,i))


#define a matrix whose columns will represent the possible basis data
data = np.transpose(np.asarray(list(map(list, itertools.product([0, 1], repeat=n)))))
data = data[:,1:len(data[0,:])]


#now check the possible combinations of basis elements to see which ones actually form a basis
if n == 2:
    for i in range(0,len(data[0,:])):
        M = np.zeros((n,n))
        M[:,0] = data[:,i]
        for j in range(0,len(data[0,:])):
            M[:,1] = data[:,j]
            rank = np.linalg.matrix_rank(M)
            if rank == n:
                numb_lin_indep_mat = numb_lin_indep_mat + 1

            if np.linalg.det(M)%2 == 0 and rank == n:
                count_epsilon += 1


if n == 3:
    for i in range(0,len(data[0,:])):
        M = np.zeros((n,n))
        M[:,0] = data[:,i]
        for j in range(0,len(data[0,:])):
            M[:,1] = data[:,j]
            for k in range(0,len(data[0,:])):
                M[:,2] = data[:,k]
                rank = np.linalg.matrix_rank(M)

                if rank == n:
                    numb_lin_indep_mat = numb_lin_indep_mat + 1
                
                if np.linalg.det(M)%2 == 0 and rank == n:
                    count_epsilon += 1


if n == 4:
    for i in range(0,len(data[0,:])):
        M = np.zeros((n,n))
        M[:,0] = data[:,i]
        for j in range(0,len(data[0,:])):
            M[:,1] = data[:,j]
            for k in range(0,len(data[0,:])):
                M[:,2] = data[:,k]
                for l in range(0,len(data[0,:])):
                    M[:,3] = data[:,l]
                    rank = np.linalg.matrix_rank(M)

                    if rank == n:
                        numb_lin_indep_mat = numb_lin_indep_mat + 1

                    if np.linalg.det(M)%2 == 0 and rank == n:
                        count_epsilon += 1



if n == 5:
    for i in range(0,len(data[0,:])):
        M = np.zeros((n,n))
        M[:,0] = data[:,i]
        for j in range(0,len(data[0,:])):
            M[:,1] = data[:,j]
            for k in range(0,len(data[0,:])):
                M[:,2] = data[:,k]
                for l in range(0,len(data[0,:])):
                    M[:,3] = data[:,l]
                    for m in range(0,len(data[0,:])):
                        M[:,4] = data[:,m]
                        rank = np.linalg.matrix_rank(M)

                        if rank == n:
                            numb_lin_indep_mat = numb_lin_indep_mat + 1

                        if int(np.linalg.det(M))%2 == 0 and rank == n:
                            count_epsilon += 1


#avoid overcounting due to counting permutations of the same basis as distinct bases
N = numb_lin_indep_mat/math.factorial(n)

#display results
print('The columns of the matrix \n',data,'\nare the possible basis elements for R ^ ',n,'consisting of only 0\'s and 1\'s.\n')
print('The number of bases that can be formed for R ^',n,'using only vectors consisting of 0\'s and 1\'s is',N)
print('\n|GL_n(F_2)|/n! for n =',n,'is',order_GL/math.factorial(n))
#print('\nThe actual \'epsilon\' term for the n =',n,'case is ',N - order_GL/math.factorial(n))
print('\nThe computed epsilon term for the n =',n,'case is ',count_epsilon/math.factorial(n))



