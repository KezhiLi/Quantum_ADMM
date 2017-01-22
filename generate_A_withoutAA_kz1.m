function [A]=generate_A_withoutAA_kz1(eta,N,P,qubit,M)

rows=randperm(N*P);
select_row=sort(rows(1:M))-1;%将d^2随机打乱，再取出前m个并排序
clear rows;
A0=[1,0;0,1];
A1=[0,1;1,0];
A2=[0,-1i;1i,0];
A3=[1,0;0,-1];
A4=[A0;A1;A2;A3];


AA=cell(1,M);

for i=1:M
    i

    str=dec2bin(select_row(i),2*qubit);
    str_4=zeros(1,qubit);
    
    Kronecker=1;
    for j=1:qubit
        str_4(j)=str2num(str(2*j))+2*str2num(str(2*j-1));
        Kronecker=kron(Kronecker,A4(2*str_4(j)+1:2*str_4(j)+2,:));
    end

    AA{i}=sparse(reshape(Kronecker,1,N*P));

    clear temp;
    clear Kronecker;

end
% if mod(i,500)~=0
%    A = sparse([A;cat(1, AA{1:mod(i,blk)})]); 
% end
size(AA)

%save(filename,variables,'-append')

A = sparse(cat(1, AA{:}))/(N^0.5);
clear AA{:}
%A = sparse((reshape(sparse([AA{1:M}]),N*P,M))');

end