function [M,outlier,X_true]=generate_rho_outlier(N,P,R,eta,nu)
LL=randn(N,R)+1i*randn(N,R);
RR=LL';
X_true=LL*RR;
while rank(X_true)<R;
LL=randn(N,R)+1i*randn(N,R);
RR=LL';
X_true=LL*RR;
end
X_true=X_true/trace(X_true);
M=ceil(eta*N*P);

outlier=randn(M,1); %M�У�1��
outlier=outlier/std(outlier); 
outlier=outlier-mean(outlier); 
a=0; 
b=sqrt(nu*norm(X_true)); 
outlier=a+b*outlier; %������ֵΪ0������Ϊ5�ĸ�˹������ 
end