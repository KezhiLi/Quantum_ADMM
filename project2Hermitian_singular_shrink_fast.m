function [ rho ] = project2Hermitian_singular_shrink_fast(rho ,alp,d,svd_sNo)
rho=reshape(rho,d,d);
% Hermitian matrix
rho=(rho+rho')/2;

% fast svd
[V,D,S]=rsvd(rho,svd_sNo);

% shrinkage
diagD=diag(D);
diagDD=max(diagD-alp,0);
diagD=diagDD+min(diagD+alp,0);

rho=V*diag(diagD)*S';

% trace of density matrix is 1
if abs(trace( rho))>1e-6
    rho = rho/((abs(trace(rho))*0.9+1*0.1));
    %rho = rho/(abs(trace(rho)));
end
%
rho=reshape(rho,d*d,1);

end

