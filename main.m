% codes for paper
% "Fast Reconstruction of High-qubit Quantum States via Low Rate Measurements"
% 
% 
% 
% 
% Copyrights to Dr. Kezhi Li, @Imperial College, 22nd Jan. 2017

clear all
clc

qubit=9;
R=1;
%SNR=40;

switch R
    case 1
        disp('R = 1')
        svd_sNo = 200;      %  adaptive parameter, keep how many singular values
        shrink_para = 1;   % is the no. of singular value decreasing
    case 2
        disp('R = 2')
        svd_sNo = 150;       
        shrink_para = 1;   
end

N=2^qubit;
P=2^qubit;

switch qubit
    case 8
        disp('8')
        load A_8.mat;
        eta = 0.03;  % 3%
    case 9
        disp('9')
        load A_9.mat;
        eta = 0.017;  % 1.7%
    case 10
        disp('10')
        load A_10.mat;
        eta = 9438/(N*P);  % 1%
    case 11
        disp('11')
        load A_11.mat;
        eta = 25166/(N*P);  % 0.6%
    case 12 
        disp('12')
        load A_12.mat;
        eta = 62500/(N*P);  % 0.3%
end

t=0.9; % parameter
gamma=1e-4; 
lamda=15;
c=1.099;%ADMM parameter
nu=0;
tol_1=1e-6; 

% generate X_true matrix, and the number of measurements M
[M,outlier,X_true]=generate_rho_outlier(N,P,R,eta,nu);

% This step generates the sensing matrix 'A'. If 'A' has been generated, we
% can load it directly.
%[A]=generate_A_withoutAA_kz1(eta,N,P,qubit,M);

% the magnitudes of noises
sigma=1e-4*norm(X_true,'fro');

% generate the measurements after corrupted by noises
b= A*(reshape(X_true,N*P,1))+sigma*randn(M,1);

% maximum iteration no.
maxite=100;
tic
% main function to reconstruct the density matrix rho
[ rho1,result_min,result_rho,result_resid,t_095,fide_095] =Robust_Quantum_fixedpoint_fast(b,A,maxite,tol_1 ,X_true,gamma,lamda,c,M,t,svd_sNo,shrink_para);  % ,svd_sNo,shrink_para
T=toc;
     
% draw results for the correspoding qubits        
xx=1:1:maxite;
figure,
plot(xx,result_rho(1,:),'r*-');
hold on
plot(xx,result_resid(1,:),'b*-');
xlabel('iteration');
ylabel('error')
title('n=');

% display the time and fidelity when 95% accuracy has been achieved
disp(t_095)
disp(fide_095)




