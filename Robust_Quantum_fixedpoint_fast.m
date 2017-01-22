 function [ rho1,result_min,result_rho,result_resid, t_095,fide_095] =Robust_Quantum_fixedpoint_fast (b,A,maxite,tol_1 ,X_true,gamma,lamda,c,M,t,svd_sNo,shrink_para)
[A_row,A_col]  = size(A);  
d = round(sqrt(A_col)); %  the nearest integers 

Li=1;

alp=t/lamda;%奇异值收缩参数
y=zeros(M,1); %y的初始值,向量
rho=zeros(d*d,1);
ite=0;
converged=0;
result_min=zeros(1,maxite);
result_rho=zeros(1,maxite);
result_resid=zeros(1,maxite);

t_095 = -1;
t_start = tic;
flag =0;
fide_095=-1;

% svd_sNo = 150;
% shrink_para = 0.95;
while ~converged
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
%     %%
%     ind_loc = abs(rho)<0.1*mean(abs(rho));
%     rho(ind_loc)=0;
%     %%
    rho1=reshape(rho,d,d);
    result_1=norm(X_true-rho1,'fro')^2/norm(X_true,'fro')^2;
    if result_1<0.05&&flag == 0
        t_095 = toc
        t_095 = [ite,t_095];
        flag = 1;
        %fide_095 = sqrt(trace(rho1*X_true)+2*sqrt(det(rho1)*det(X_true)));
        fide_095 = Fidelity(X_true,rho1); %trace(sqrt(sqrt(X_true)*rho1*sqrt(X_true)));
    end
    svd_sNo = round(shrink_para*svd_sNo);
 %   [V,D,U]=rsvd(rho1,svd_sNo);  %svdecon
    %[V,D,U]=svd(rho1);
 %   diagD=diag(D);
 %   result_2=norm(diagD,1); 
   
    result_rho(ite)=result_1;
 %   result_min(ite)=result_2;
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

