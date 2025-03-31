% Author: Maryam Bagherian
% Date: 03-31-2025
% 
% References:
% [1] Bagherian, M, et al., "Coupled tensor-tensor completion methods"
% [2] Bagherian, M, et al., "Coupled matrix–matrix and coupled tensor–matrix completion methods for predicting drug–target interactions." 
%     Briefings in bioinformatics 22, no. 2 (2021): 2161-2171, doi:
%     10.1093/bib/bbaa025.
% [3] Bagherian, M. "Tensor denoising via dual Schatten norms." Optimization Letters 18, no. 5 (2024): 1285-1301.
% This code also uses some parts of the follwoing refrence: 
% [4] Chen, Yi-Lei, Chiou-Ting Hsu, and Hong-Yuan Mark Liao. "Simultaneous tensor decomposition and completion using factor priors." 
%     IEEE transactions on pattern analysis and machine intelligence 36, no. 3 (2013): 577-591.




%% Para
% scale and epsilon manually set
function out = CTTC(X, Side, iter,m_rate)
        %% Description
        % VSet: the subsets encoding sub-manifolds, e.g., {[1,2],[3]} means that
        % the 1st and 2nd sub-manifolds are encoded simultaneously and the
        % 3rd sub-manifold is encoded individually
        % Rate: the downsampling rates of different Laplacian graphs
        % Affinity: a cell array, where the i-th element denotes an affinity matrix
        % to encode the intra-factor relation w.r.t. the i-th tensor's
        % dimension
        % X: Initial Tensor
        % Side: A Cell consiting three tensors of side informations
        % scale: vector of 3 scaling factors
        % epsilon: vector of 3 regularization paramter
        % iter: number of iterrations
        % m_rate: missing rate 
        %% 
rng(42);
iters = [50];
miss_rate = [0.1, 0.2, 0.3];
% for gap = 1:10
for m_rate = miss_rate
    for maxitr = iters 
        rse_list = [];
        g_time_list = [];
        time_list = [];
%         for gap = 1:10
%  Load the sample data
        X1 = load(".../filtered_target_tensor_1.mat"); 
        X2 = load(".../filtered_target_tensor_2.mat");
        X3 = load(".../filtered_target_tensor_3.mat");
        X4 = load(".../filtered_target_tensor_4.mat");
        X5 = load(".../filtered_target_tensor_5.mat");
        X1 = X1.labpcexport;
        X2 = X2.labpcexport;
        X3 = X3.labpcexport;
        X4 = X4.labpcexport;
        X5 = X5.labpcexport;
        X = cat(2, X1, X2, X3, X4, X5);
        ds_1 = load(".../drug_sim_first_half.mat");
        ds_2 = load(".../drug_sim_sec_half.mat");
        ds_3 = load("drug_sim_third.mat");
        gs = load("gene_sim.mat");
        cs = load("cell_sim.mat");
        ds_1 = ds_1.X;
        ds_2 = ds_2.X;
        ds_3 = ds_3.X;
        gs = gs.X;
        cs = cs.X;
        ds = cat(1, ds_1, ds_2, ds_3);
        ds = ds([1, 5], :, :);
        gs = (gs + 1)/2;
        ds = permute(ds, [2, 3, 1]);
        gs = permute(gs, [2, 3, 1]);
        cs = permute(cs, [2, 3, 1]);

        T = {gs; ds; cs};
        % m_rate = 0.2;
        X = (X + (abs(min(X(:)))))/(max(X(:)) + abs(min(X(:))));
        X(isnan(X)) = -1;
        Data=X;
        kappa=1.5849;
        omega=3.9811;tau=0.01;gnns=2; pnns=100;
        % maxitr=10; % Normally, we take it as 50
        mode_dim=1;% 0 does NOT ignore the last factor, 1 ignores the last factor
        mode_pm=0;
    
        Rate=[1,1];
        tsize=size(Data);
        % first and second factors are assessed individually over the manifold
        VSet={[1],[2]};
        Affinity={zeros(tsize(1)), zeros(tsize(2)),zeros(tsize(3))};
        para = initial_para(kappa,omega,tau,gnns,pnns,maxitr,mode_dim,...
        mode_pm,mode_noise,VSet,Rate,Affinity,tsize);
        
        
        tic;
        tsize = size(X);
        Xg = X; % copying for ground truth value
        %%%%____Creating a tensor with m_rate% missing entries___%%%%%
        [m1, n1, p1] = size(X);
        random_idx = randsample(m1, int32(m_rate * m1));
        mark = zeros(tsize);
        mark(random_idx, :, :) = 1;
        mark = logical(mark);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Initializing the side tensors 
        X1=X(:,2,3);sX1=size(X1);
        X2=X(1,:,3)';sX2=size(X2);
        X3=permute(X(1,2,:),[3,1,2]);sX3=size(X3);
        X1(1) = 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i=1:sX1(1)
                 for j=1:sX1(1)
                     Side1(i,j)=dist(X1(i),X1(j));
                 end
            end
            for i=1:sX2(1)
                 for j=1:sX2(1)
                     Side2(i,j)=dist(X2(i),X2(j));
                 end
            end
            for i=1:sX3(1)
                 for j=1:sX3(1)
                     Side3(i,j)=dist(X3(i),X3(j));
                 end
            end
        VXX={Side1,Side2,Side3};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % initialize manifold graphs
        if ~isfield(para,'VSet')
            para.H = [];
        else
            para.H = construct_graphL(tsize,para.VSet,para.Rate,para.gnns,para.Affinity);
            for i = 1 : size(para.H,1)
                if size(para.H,2)==1 || numel(para.H{i,2})==0
                    for j = 1 : numel(tsize)
                        para.Ds{i,j} = eye(tsize(j));
                    end
                else
                    for j = 1 : numel(tsize)
                        A = randn(tsize(j)); 
                        B = imresize(A,[tsize(j) round(tsize(j)*para.Rate(i))],'bilinear');
                        para.Ds{i,j} = A\B;
                    end
                end
            end
        end
        
          
        V={};
        vsize = tsize;
        N = numel(tsize);
        if para.mode_dim, N = N-1; end
        if para.mode_PoM
            randn('seed',1);
            for i = 1 : N
                V{i} = randn(tsize(i));
                V{i} = V{i}/norm(V{i});
                vsize(i) = size(V{i},1);
            end
            Xt = HaLRTC(X,mark,ones(N,1),10^-2,1.1,100,N,Xg);

            P = initial_PoM(Xt,para.gnns^2,para.pnns,N);
        else
        
            for i = 1 : N
                V{i} = eye(tsize(i));
               % Distanc Metric Learning Loop using CMMC function
               % V{i}=CMMC(VXX, V{i}, VYY, xscale, yscale, xepsilon, 
               %      ... yepsilon,num_iter);
               % INPUTS (8)
                % VXX (ixi double): The matrix to be coupled along the X-axis
                % V{i} (ixj double): The matrix to be completed by coupling MXX and MXY.
                % VYY (jxj double): The matrix to be coupled along the Y-axis.
                % xscale (double): The scaling factor for VXX
                % yscale (double): The scaling factor for VYY
                % xepsilon (double): Additional value to add in the diagnoal elements of the
                %   VXX
                % yepsilon (double): Additional value to add in the diagnoal elements of the
                %   VYY
                
                % OUTPUTS (1)
                % V{i} (nxn double): The completed input V{i}.

                % NOTE: 
                % To execute CMMC Function one needs to place the function CoupledTensorCompletion
                % in the same folder. 
                
                % Citation: Bagherian et al. (2020) Coupled Matrix--Matrix and
                % Tensor--Matrix Completion Methods for Predicting Drug--Target
                % Interactions. Briefings in Bioinformatics
              
                V{i}=CMMC(VXX{i}, V{i}, VXX{i}, 0.1, 0.1, 0.1, 0.1, 20);
                vsize(i) = size(V{i},1);
                P{i} = (1:tsize(i))';
            end
             
        end
        % initialize core tensor Z as the Mahalanobis Distance 
        Z = X;
        Y = zeros(tsize);
        % initialize algorithm parameters
        norm_gt = norm(Xg(:));
        norm_x = norm(X(:));
        para.alpha = ones(N,1);
        para.gamma = (para.omega/para.tau)/(norm_x^2);
        xxt = reshape(X,tsize(1),[]);
        xxt = norm(xxt*xxt');
        ita = 1/(para.tau*xxt);
        ct = zeros(1,N);
        for i = 1 : size(para.H,1)
            ct = ct+double(para.VSet{i});
        end
        for i = 1 : size(para.H,1)
            list = 1 : N;
            list = list(para.VSet{i});
            for j = 1 : size(para.H,2)
                U{j} = para.Ds{i,list(j)};
            end
            for j = 1 : size(para.H,2)
                llt = reshape(TensorChainProduct(para.H{i,j},U,[1:j-1 j+1:size(para.H,2)]),tsize(list(j)),[]);
                llt = norm(llt*llt');
                para.H{i,j} = para.H{i,j}*para.kappa*sqrt(ita*xxt/(2*llt*ct(list(j))));
            end
        end
        % message
        g_time_list = [g_time_list;toc];
        
        disp('------------------------------------------------------------------------------');
        disp('--                          Start CTTC algorithm..                          --');
        disp('------------------------------------------------------------------------------');
        % Main algorithm
        lambda = ita*(1.1^(para.maxitr))/2;
        tic;
        for itr = 1 : para.maxitr
            
            % update V1,...,Vn as local metrics
            [V,P,rank_vi,pstr] = optimize_otv_V(X,Y,Z,V,P,ita,tsize,vsize,para);
            % update Mahalanobis Distance
            Z = optimize_Z(V,X,Y,ita,para.gamma);
            % update the initial tensor with missing entries through the
            % distance functions
            Xt = TensorChainProductT(Z,V,1:numel(V));
            X(mark) = Xt(mark)-Y(mark)/ita;
            if para.mode_nse
                X(~mark) = ((ita*Xt(~mark)-Y(~mark))+lambda*Xg(~mark))/(ita+lambda);
            end
            residual = (norm(X(:)-Xt(:)))/norm(Xt(:));
            % update Y
            Y = Y+ita*(X-Xt);
            % assessment
            info.rse(itr) = norm(X(mark)-Xg(mark))/norm_gt;
            info.fit(itr) = 1 - norm(X(mark)-Xg(mark))/norm_gt;
            info.rank_vi(:,itr) = rank_vi;
            info.residual(:,itr) = residual;

            % stopping criterion
            if residual<=0.00001, break; end
            ita = ita*1.2;
            if itr>=50
                ita = ita*2;
            end
            
        end
        Core = Z;
        for i = 1 : N
            V{i} = V{i}';
        end
        info.P = P;
        rse_list = [rse_list; info.rse(itr)];

        fprintf("mrate: %d; maxitr: %d; rse: %d\n", m_rate, maxitr, info.rse(itr));

        elapsedTime = toc;  % Stop timer and get the elapsed time
        time_list = [time_list; elapsedTime];
        fprintf('Elapsed time: %.4f seconds\n', elapsedTime);
    end
        
        
end

function para = initial_para(kappa,omega,tau,gnns,pnns,maxitr,mode_dim,mode_PoM,mode_nse,VSet,Rate,Affinity,tsize)
        % Description
        % VSet: the subsets encoding sub-manifolds, e.g., {[1,2],[3]} means that
        % the 1st and 2nd sub-manifolds are encoded simultaneously and the
        % 3rd sub-manifold is encoded individually
        % Rate: the downsampling rates of different Laplacian graphs
        % Affinity: a cell array, where the i-th element denotes an affinity matrix
        % to encode the intra-factor relation w.r.t. the i-th tensor's
        % dimension
        % Initialization
        para.kappa = kappa; % the weight of MGE
        para.omega = omega; % (||Z||_F)^2 ~= omega*(ita*(||X-Zx1V1x2...xnVn||_F)^2)/2
        para.tau = tau; % the 1st threshold for rank minimization
        para.gnns = gnns; % # effective neighbors on manifold graphs
        para.pnns = pnns; % # possible destinations on manifold graphs
        para.maxitr = maxitr; % maximal iteration of IALM algorithm
        para.mode_dim = mode_dim; % ignore the last dimension (true) or not (false)
        para.mode_PoM = mode_PoM; % update permutation matrices (true) or not (false)
        para.mode_nse = mode_nse; % ignore the observation noise (false) or not (true)
        % convert graph indices into binary indicator
        num_g = numel(VSet);
        N = numel(tsize);
        id=[];
        if mode_dim, N=N-1; end
        for i = 1 : num_g
        id{i} = logical(zeros(1,N));
        if sum(VSet{i}>N)>0
        error('Wrong index of manifold graph!');
        end
        id{i}(VSet{i}) = logical(1);
        end
        para.VSet = id;
        % ensure # downsampling rates == # graphs
        if numel(Rate)~=num_g
        disp('# graphs is incorrect...');
        disp('Correct it automatically...');
        else
        Rate(num_g+1:end) = [];
        Rate(end+1:num_g) = 1;
        end
        para.Rate = Rate;
        % confirm the affinity matrix
        for i = 1 : N
        para.Affinity{i} = Affinity{i};
  
        end
        end
end