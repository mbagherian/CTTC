function MXY = CTMC(TXX, MXY, TYY, xscale, yscale, xepsilon, yepsilon,...
    num_iter)
% TXX=T{1}; MXY=Mi;TYY=TT{1};xscale=0.6; yscale=0.6; xepsilon=0.8;
% yepsilon=0.8; num_iter=1;


% INPUTS (8)
% TXX (ixixn double): The tensor to be coupled along the X-axis
% MXY (ixj double): The matrix to be completed by coupling MXX and MXY.
% TYY (jxjxm double): The tensor to be coupled along the Y-axis.
% xscale (double): The scaling factor for TXX
% yscale (double): The scaling factor for TYY
% xepsilon (double): Additional value to add in the diagnoal elements of the
%   TXX
% yepsilon (double): Additional value to add in the diagnoal elements of the
%   TYY

% OUTPUTS (1)
% MXY (ixj double): The completed input MXY.

% Version: MATLAB R2018b
% OS: Windows 10
% Author: Renaid Kim
% Date: 12/2/2019

% Citation: Bagherian et al. (2020) Coupled Matrix--Matrix and
% Tensor--Matrix Completion Methods for Predicting Drug--Target
% Interactions. Briefings in Bioinformatics


% Preprocess TXX/TYY

for k = 1:size(TXX, 3)
    slice = TXX(:,:,k);
    slice = slice^2;
    Diag=diag(diag(slice).^(-1/2));
    slice=Diag*slice*Diag;
    TXX(:,:,k) = slice;
end

for k = 1:size(TYY, 3)
    slice = TYY(:,:,k);
    slice = slice^2;
    Diag=diag(diag(slice).^(-1/2));
    slice=Diag*slice*Diag;
    TYY(:,:,k) = slice;
end

a = length(TXX);
b = length(TYY);

% create inputs for coupled tensor completion.
C_1 = zeros(a);
for k = 1:size(TXX, 3)
    slice = TXX(:,:,k);
    C_slice = xscale*slice+xepsilon*eye(size(slice));
    C_1 = C_1 + C_slice;
end

C_2 = zeros(b);
for k = 1:size(TYY, 3)
    slice = TYY(:,:,k);
    C_slice = yscale*slice+yepsilon*eye(size(slice));
    C_2 = C_2 + C_slice;
end

C{1} = C_1;
C{2} = C_2;

Omega=[];
v=[];
W = MXY;

for i=1:a
    for j=1:b
        if W(i,j) == 1
            Omega=[Omega;i,j];
            v=[v;W(i,j)];
        end
    end
end

[s,~]=size(Omega);
%Diag={};
[Diag,u]=CoupledTensorCompletion(Omega,v,C,num_iter);

E=zeros(a,b);

for i=1:s
    E(Omega(i,1),Omega(i,2))=u(i);
end

MXY=double(Diag{1})*E*double(Diag{2});
end

