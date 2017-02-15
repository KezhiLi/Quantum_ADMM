function [M,outlier,X_true]=generate_rho_outlier(N,P,R,eta,nu)
% function that generates true density matrix, outliers, and number of
% measurements M
% 
% 
% Copyrights to Dr. Kezhi Li, @Imperial College, 22nd Jan. 2017

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

outlier=randn(M,1); % generate random outliers 
outlier=outlier/std(outlier); 
outlier=outlier-mean(outlier); 
a=0; 
b=sqrt(nu*norm(X_true)); 
outlier=a+b*outlier; 
end