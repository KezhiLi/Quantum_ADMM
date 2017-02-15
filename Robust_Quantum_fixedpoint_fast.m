function [ rho1,result_min,result_rho,result_resid, t_095,fide_095] = Robust_Quantum_fixedpoint_fast (b,A,maxite,tol_1 ,X_true,gamma,lamda,c,M,t,svd_sNo,shrink_para)
% function that generates the reconstruction result
% rho1: the reconstructed density matrix rho
% result_rho: reconstructed accuracy
% t_095: the time uses when the algorithm acheives 95% accuracy
% fide_095: the time uses when the algorithm acheives 95% fidelity
%
% Copyrights to Dr. Kezhi Li, @Imperial College, 22nd Jan. 2017

[A_row,A_col]  = size(A);
d = round(sqrt(A_col)); %  the nearest integers

alp=t/lamda;
y=zeros(M,1);
rho=zeros(d*d,1);
ite=0;
converged=0;
result_min=zeros(1,maxite);
result_rho=zeros(1,maxite);
result_resid=zeros(1,maxite);

t_095 = -1;
flag =0;
fide_095=-1;
tic;

while ~converged
    % display iteration
    ite=ite+1
    
    if mod(ite,10)==0
        disp(['iterations:' num2str(ite)]);
    end    
    xi=(gamma*lamda/(1+gamma*lamda))*(-y/lamda-(A*rho-b));
    
    rho=rho-t*(A'*(A*rho+xi-b+y/lamda));%向量
    rho=project2Hermitian_singular_shrink_fast(rho,alp,d,svd_sNo); %向量
    
    resid=A*rho+xi-b;
    y=y+c*lamda*resid;
    
    result_3=norm(resid,2);

    rho1=reshape(rho,d,d);
    result_1=norm(X_true-rho1,'fro')^2/norm(X_true,'fro')^2;
    if result_1<0.05&&flag == 0
        t_095 = toc
        t_095 = [ite,t_095];
        flag = 1;
        fide_095 = Fidelity(X_true,rho1); %trace(sqrt(sqrt(X_true)*rho1*sqrt(X_true)));
    end
    svd_sNo = round(shrink_para*svd_sNo);
    
    result_rho(ite)=result_1;
    result_resid(ite)=result_3;
    
    stop=norm(resid,'fro');
    
    if stop<tol_1
        converged=1;
    end
    if ~converged && ite>=maxite
        disp('maximum iteration reached');
        converged=1;
    end
end
end

