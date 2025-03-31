function MXY = CMMC(MXX, MXY, MYY, xscale, yscale, xepsilon, yepsilon,...
    num_iter)
% INPUTS (8)
% MXX (ixi double): The matrix to be coupled along the X-axis
% MXY (ixj double): The matrix to be completed by coupling MXX and MXY.
% MYY (jxj double): The matrix to be coupled along the Y-axis.
% xscale (double): The scaling factor for MXX
% yscale (double): The scaling factor for MYY
% xepsilon (double): Additional value to add in the diagnoal elements of the
%   MXX
% yepsilon (double): Additional value to add in the diagnoal elements of the
%   MYY

% OUTPUTS (1)
% MXY (nxn double): The completed input MXY.

% Version: MATLAB R2018b
% OS: Windows 10
% Author: Renaid Kim
% Date: 11/26/2019

% Citation: Bagherian et al. (2020) Coupled Matrix--Matrix and
% Tensor--Matrix Completion Methods for Predicting Drug--Target
% Interactions. Briefings in Bioinformatics


% Preprocess MXX/MXY
MXX=MXX^2;
D=diag(diag(MXX).^(-1/2));
MXX=D*MXX*D;
MYY=MYY^2;
D=diag(diag(MYY).^(-1/2));
MYY=D*MYY*D;

a = length(MXX);
b = length(MYY);

% create inputs for coupled tensor completion.
C{1}=xscale*MXX+xepsilon*eye(size(MXX));
C{2}=yscale*MYY+yepsilon*eye(size(MYY));

Omega=[];
v=[];
W = MXY;
for i=1:a
    for j=1:b
        if W(i,j)==1
            Omega=[Omega;i,j];
            v=[v;W(i,j)];
        end
    end
end

[s,~]=size(Omega);

[D,u]=CoupledTensorCompletion(Omega,v,C,num_iter);

E=zeros(a,b);

for i=1:s
    E(Omega(i,1),Omega(i,2))=u(i);
end

MXY=D{1}*E*D{2};
end