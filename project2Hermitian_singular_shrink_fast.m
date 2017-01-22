function [ rho ] = project2Hermitian_singular_shrink_fast( rho ,alp,d,svd_sNo)
rho=reshape(rho,d,d);
rho=(rho+rho')/2;
 
%[V,D,S]=svd(rho);
 [V,D,S]=rsvd(rho,svd_sNo);
 %diagD=diag((max(abs(D)-tao,0)).*sign(D));
 
 
 diagD=diag(D);
 diagDD=max(diagD-alp,0);
 diagD=diagDD+min(diagD+alp,0);

 rho=V*diag(diagD)*S';
  %
 if abs(trace( rho))>0.000001
    rho = rho/((abs(trace(rho))*0.9+1*0.1));
    %rho = rho/(abs(trace(rho)));
 end
 %
 rho=reshape(rho,d*d,1);
 
end

