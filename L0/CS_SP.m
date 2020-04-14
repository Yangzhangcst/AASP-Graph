function [A]=CS_SP(D,X,L)
%=============================================
% Sparse coding of a group of signals based on a given 
% dictionary and specified number of atoms to use. 
% input arguments: 
%       D - the dictionary (its columns MUST be normalized).
%       X - the signals to represent
%       L - the max. number of coefficients for each signal.
% output arguments: 
%       A - sparse coefficient matrix.
%=============================================
[n,P]=size(X);
[n,K]=size(D);
for k=1:1:P,
    a=[];
    x=X(:,k);
    residual=x;
    indx=zeros(L,1);
    for j=1:1:L
        proj=D'*residual;
        [~,pos]=sort(abs(proj),'descend');
        indx = pos(1:j);
        a=pinv(D(:,indx))*x;
        [~,pos]=sort(abs(a),'descend');
        indx = indx(pos(1:j));
        a = a(pos(1:j));
        residual=x-D(:,indx)*a;
        if sum(residual.^2)<1e-6,  break;   end
    end;
    temp=zeros(K,1);
    temp(indx)=a;
    A(:,k)=sparse(temp);
end;
return;