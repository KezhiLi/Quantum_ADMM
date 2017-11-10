# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 15:01:38 2017

@author: kezhili
"""
import scipy.io as sio
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import math, time, sklearn.decomposition
from scipy import sparse

qubit = 8
R = 1
N=2**qubit;
P=2**qubit;

if R == 1:
    svd_sNo = 150
    shrink_para = 1
elif R == 2:
    shrink_para = 1 
    
'''
function [M,outlier,X_true]=generate_rho_outlier(N,P,R,eta,nu)
'''
def generate_rho_outlier(N,P,R,eta,nu):
# function that generates true density matrix, outliers, and number of measurements M
# 
# Copyrights to Dr. Kezhi Li, @Imperial College, 22nd Jan. 2017

    LL=np.random.randn(N,R)+1j*np.random.randn(N,R)
    RR=LL.T
    X_true=np.dot(LL,RR)
    while LA.matrix_rank(X_true)<R:
        LL=np.random.randn(N,R)+1j*np.random.randn(N,R)
        RR=LL.T
        X_true=np.dot(LL,RR)

    X_true=X_true/np.trace(X_true)
    M=math.ceil(eta*N*P)
    
    outlier=np.random.randn(M,1) # generate random outliers 
    outlier=outlier/np.std(outlier) 
    outlier=outlier-np.mean(outlier)
    a=0
    b=np.sqrt(nu*LA.norm(X_true))
    outlier=a+b*outlier

    return M,outlier,X_true

'''
function Fidelity
'''
def Fidelity(rho,sigma):
    sq_rho,res = LA.sqrtm(rho); # need "res" parameter to suppress MATLAB singularity warning
    sq_fid,res = LA.sqrtm(np.dot(np.dot(sq_rho,sigma),sq_rho));
    fid = np.real(np.trace(sq_fid)); # finally, compute the fidelity
    return fid
'''
function [ rho1,result_min,result_rho,result_resid, t_095,fide_095] = Robust_Quantum_fixedpoint_fast (b,A,maxite,tol_1 ,X_true,gamma,lamda,c,M,t,svd_sNo,shrink_para)
'''
def Robust_Quantum_fixedpoint_fast_sq(b,A,maxite,tol_1 ,X_true,gamma,lamda,c,M,t,svd_sNo,shrink_para):
# function that generates the reconstruction result
# rho1: the reconstructed density matrix rho
# result_rho: reconstructed accuracy
# t_095: the time uses when the algorithm acheives 95% accuracy
# fide_095: the time uses when the algorithm acheives 95% fidelity
#
# Copyrights to Dr. Kezhi Li, @Imperial College, 22nd Jan. 2017

    A_col  = A.shape[1]*A.shape[2]
    d = int(round(np.sqrt(A_col))) #  the nearest integers
    
    alp=t/lamda
    y=np.zeros(M)
    rho_sq=np.zeros((d,d))
    ite = 0
    converged = 0
    result_min=np.zeros(maxite)
    result_rho=np.zeros(maxite)
    result_resid=np.zeros(maxite)
    
    t_095 = list([-1])
    flag =0
    fide_095=-1
    tic = time.time()

    while ~converged:
        # display iteration
        ite=ite+1
        print(ite)
        
        if ite%10==0:
            print("iterations:"+str(ite))
        
        AA_rho0 = np.einsum('ijk,jk',AA,rho_sq)
        #xi=(gamma*lamda/(1+gamma*lamda))*(-y/lamda-(np.dot(A,rho0)-b))    
        xi=(gamma*lamda/(1+gamma*lamda))*(-y/lamda-(AA_rho0-b))    
        #rho1=rho0-t*np.dot(A.T,(np.dot(A,rho0)+xi-b+y/lamda))
        rho1_sq=rho_sq-t*np.einsum('ijk,i',AA,(AA_rho0+xi-b+y/lamda))
        
        rho_sq=project2Hermitian_singular_shrink_fast_sq(rho1_sq,alp,d,svd_sNo)
        
        #resid=np.dot(A,rho)+xi-b
        resid=np.einsum('ijk,jk',AA,rho_sq)+xi-b
        y=y+c*lamda*resid
        
        result_3=LA.norm(resid,2)
    
        result_1=(LA.norm(X_true-rho_sq,'fro'))**2/(LA.norm(X_true,'fro'))**2
        if result_1<0.05 and flag == 0:
            t_095 = [time.time() - tic]
            t_095.insert(0,ite)
            flag = 1
            fide_095 = Fidelity(X_true,rho_sq); #trace(sqrt(sqrt(X_true)*rho1*sqrt(X_true)))

        svd_sNo = round(shrink_para*svd_sNo)
        
        result_rho[ite-1]=result_1
        result_resid[ite-1]=result_3
        
        stop=LA.norm(resid)
        
        if stop<tol_1:
            converged=1

        if ~converged and ite>=maxite:
            print('maximum iteration reached')
            converged=1

    return rho1_sq,result_min,result_rho,result_resid, t_095,fide_095

'''
main start    
'''
def project2Hermitian_singular_shrink_fast_sq(rho_sq,alp,d,svd_sNo):
    # function that carries out the singular shrink step 
    # 
    # 
    #%   Copyrights to Dr. Kezhi Li, @Imperial College, 22nd Jan. 2017
    
    
    # Hermitian matrix
    rho=(rho_sq+rho_sq.T)/2
    
    # fast svd
    tsvd = sklearn.decomposition.TruncatedSVD(svd_sNo, algorithm="randomized", n_iter=1)
    U, D, VT=tsvd.fit_transform2(rho)
    
    # shrinkage
#    diagDD = list(map(max, zip(diagD-alp, np.zeros(len(diagD)))))
#    diagD = diagDD + list(map(min, zip(diagDD+alp, np.zeros(len(diagDD)))))
    diagDD = np.fmax(D-alp, np.zeros(len(D)))
    diagD = diagDD + np.fmax(diagDD+alp, np.zeros(len(diagDD)))
    
    rho=np.dot(np.dot(U,np.diag(diagD)),VT)
    
    # trace of density matrix is 1
    if abs(np.trace(rho))>1e-6:
        rho = rho/((abs(np.trace(rho))*0.9+1*0.1))
        #rho = rho/(abs(trace(rho)))
    return rho
'''
main start    
'''
num = 8

if num == 8:
    print(8)
    data = sio.loadmat('A_8.mat')
    A = data['A']
    eta = 0.03; 
elif num == 9:
    print(9)
    data = sio.loadmat('A_9.mat')
    A = data['A']
    eta = 0.017; 
elif num == 10:
    print(10)
    data = sio.loadmat('A_10.mat')
    A = data['A']
    eta = 9438/(N*P); 
elif num == 11:
    print(11)
    data = sio.loadmat('A_11.mat')
    A = data['A']
    eta = 25166/(N*P); 
elif num == 12:
    print(12)
    data = sio.loadmat('A_12.mat')
    A = data['A']
    eta = 62500/(N*P); 

t=0.9
gamma=1e-4
lamda=15
c=1.099
nu=0
tol_1=1e-6

AA = np.reshape(A.toarray(),(A.shape[0],N,P))
# generate X_true matrix, and the number of measurements M
M,outlier,X_true=generate_rho_outlier(N,P,R,eta,nu)

# This step generates the sensing matrix 'A'. If 'A' has been generated, we
# can load it directly.
#[A]=generate_A_withoutAA_kz1(eta,N,P,qubit,M);

# the magnitudes of noises
sigma=1e-4*LA.norm(X_true,'fro')

# generate the measurements after corrupted by noises
#b= np.dot(A,(np.reshape(X_true,(N*P,1)))) # +sigma*np.random.randn(M)
b = np.einsum('ijk,jk',AA,X_true)

# maximum iteration no.
maxite=100
# main function to reconstruct the density matrix rho
rho1_sq,result_min,result_rho,result_resid,t_095,fide_095 =Robust_Quantum_fixedpoint_fast_sq(b,AA,maxite,tol_1 ,X_true,gamma,lamda,c,M,t,svd_sNo,shrink_para) # ,svd_sNo,shrink_para
     
# draw results for the correspoding qubits        
xx=np.array(range(maxite))
plt.plot(xx,result_rho[0,:],'r--')
plt.plot(xx,result_resid[0,:],'b^')
plt.xlabel('iteration')
plt.ylabel('error')
plt.show()

# display the time and fidelity when 95% accuracy has been achieved
print(t_095)
print(fide_095)


