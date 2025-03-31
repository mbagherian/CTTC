function [B,u] = CoupledTensorCompletion(Omega,v,C,iter )
% Omega: an m x d matrix, where m is the number of known entries, and d is
% the number of modes in the tensor
%
% v: an m x 1 vector of prescribed values for the tensor
% 
% C: m x 1 cell where C{i} is a regularization matrix for mode i
%
% iter: number of iterations
B={}; 
[m,d]=size(Omega);
for i=1:d
    n(i)=size(C{i},1);
end


% we start initialization
for i=1:d
    B{i}=eye(n(i));
    Z(:,:,i)=B{i}(Omega(:,i),Omega(:,i));
end
u=v;

% we start the main loop
for l=1:iter
 %   l
    for i=1:d
        Y=prod(Z(:,:,[1:i-1,i+1:d]),3);
 
    % we will update A{i} for i=1,...,d
        M=zeros(n(i),n(i));
        for j=1:m
            for k=1:m
                x=Omega(j,i);
                y=Omega(k,i);
                M(x,y)=M(x,y)+u(j)*u(k)*Y(j,k);
            end
        end
        % update A{i}
        B{i}=(C{i}+B{i}.*M.*B{i}');
        B{i}=double(B{i});
        dt=rcond(B{i})^(-1/n(i));
        B{i}=dt*B{i};
        B{i}=tensor(B{i});
        Z(:,:,i)=B{i}(Omega(:,i),Omega(:,i));
        X=Y.*Z(:,:,i);
        u=X\v; % update u
    end  
end

