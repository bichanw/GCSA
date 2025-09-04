function [U_vi_cell,V_vi_cell,to_save] = multiRRR_wARD(Xin,Yout,rnks,lam0,opts,ard_params)
% [wU,wVt,wwfilts,fval] = bilinearMultifiltRRR_coordAscent(Xin,Yout,rnks,lambda,opts)
% 
% Computes low-rank regression estimate using coordinate ascent in the
% multi-filter reduced rank regression (RRR) problem setting 
% the parameter vector.
%
% Finds solution to: 
%
%   argmin_{Wi} ||Y - \sum_i Xi Wi||_F^2 + lambda*\sum_i ||W_i||_F^2
%
% where the matrices {Wi} are all low rank, parametrized as
%
%   Wi = Ui Vi^T, 
%
% and Ui and Vi^T denote row and column vector matrices, respectively.
%
% Inputs:
% -------
%    Xin [{1 x k}]  = cell array of input population responses {(T x nin1), ..., (T x nink)}
%   Yout [T x nout] = output population response matrix 
%   rnks [1 x k]    = rank of each low-rank filter
%
%   lambda [1 x 1]  = ridge parameter (optional)  
%
%     opts [struct] =   options struct (optional)
%         fields: 'MaxIter' [25], 'TolFun' [1e-6], 'Display' ['iter'|'off']
%
% Outputs:
% -------
%        wU [{k x 1}] = cell array of column vector matrices 
%       wVt [{k x 1}] = cell array of row vector matrices    
%   wwfilts [{k x 1}] = cell array of low rank filters
%      fval [1 x 1]   = loss function (squared error)


% Update V^T and compute the covariance over V^T given the current U (which has been set to a semi-orthogonal matrix using SVD).
% Compute one update for theta and sig^2 using that posterior mean and covariance.
% Convert V^T to semi-orthogonal using SVD. 
% Update U using fixed V^T (which is semi-orthogonal).  Note that the ridge penalty is now lambda = theta*sigma^2.  (Should check to confirm that is correct). 
% Set U to semi-orthogonal using SVD


% initialization for function testing
% Xin  = {squeeze(spk{1}(end,:,:))};
% Yout = squeeze(spk{2}(end,:,:));
% rnks = 4;
% lam0 = 1e3;
% opts = struct('MaxIter',1e3,'TolFun',1e-6,'Display','iter');

% code that generate simulation data, deactivated 
if false
    switch 1
    case 1
        % rrr simu
        [X_fit,Yout,U,V,params] = simu.RRR(struct('signse',1,'magnitude',0.7,'nx',10,'rnk',3));
        flts = struct('wfilt',U*V');
        Xin = {X_fit};
        % [rnk1.U,rnk1.V,rnk1_hprs] = rankr_VI({X_fit},Yout,params.rnk,0,struct('MaxIter',1e4,'TolFun',1e-6));
        % plt_cmp_matrices({rnk1.U*rnk1.V,U_vi_cell{1}*V_vi_cell{1}});ef;
    case 2
        % ar1 simu
        [X_fit,Y_hist_fit,Yout,tr_ind,flts,params] = simu.AR1(struct('rnk',[2 3],'signse',.3));
        Xin = {X_fit,Y_hist_fit{1}};
    end
    rnks = params.rnk;
    lam0 = 0;
    opts = struct('MaxIter',1e3,'TolFun',1e-6,'Display','iter','Init','svd');
    ard_params = params2ard_true(flts,rnks);
end


% if we're fitting the parameters or using true ones
if_fit_U = true;
if_fit_V = true;
if_fit_theta = true;


% ---------------------------------------------------
% set optimization options
% ---------------------------------------------------
if (nargin < 5) || isempty(opts)
    opts.default = true;
end
if ~isfield(opts, 'MaxIter'); opts.MaxIter = 25; end
if ~isfield(opts, 'TolFun'); opts.TolFun = 1e-6; end
if ~isfield(opts, 'Display'); opts.Display = 'iter'; end
if ~isfield(opts, 'Init'); opts.Init = 'svd'; end


% ---------------------------------------------------
% Extract sizes of inputs
% ---------------------------------------------------

nin = cellfun(@(x)size(x,2),Xin); % get # of cells in each population
ninpops = length(nin); % number of input populations
nintot = sum(nin);     % total number of input neurons
nout = size(Yout,2);
nsamp = size(Yout,1);
nU = sum(nin.*rnks);
nV = sum(rnks)*nout;

% check inputs
if length(rnks)~=ninpops
    error('length of ``rnks'' doesn''t match # of input populations (# cells in Xin)');
end

% ---------------------------------------------------
% Compute Sufficient Statistitics
% ---------------------------------------------------

% Compute sufficient statistics
Xfull = cell2mat(Xin);
XX = (Xfull'*Xfull);
% if length(lambda) == 1 % uniform regularization strength
%      + lambda*speye(nintot); 
% elseif length(lambda) == ninpops % different regularization strengths for each input population
%     XX = (Xfull'*Xfull) + sparse(1:nintot,1:nintot,cell2mat(arrayfun(@(ii) lambda(ii)*ones(1,nin(ii)),1:ninpops,'uni',0))); 
% else
%     error('lambda must be a scalar or a vector of length equal to the number of input populations');
% end
XY = Xfull'*Yout; % input-output cross-covariance

% Divide sufficient statistics into blocks of a cell array
% XXc = mat2cell(XX,nin,nin);    % divide XX into blocks
% XXc_cols = mat2cell(XX,nintot,nin); % divide XX into columns 
% XYc = mat2cell(XY,nin,nout);   % divide XY into blocks
% rnks_cell = num2cell(rnks(:)); % filter ranks as cell array
totrnks = sum(rnks);
prod_rnk_nin = rnks .* nin; 

% ---------------------------------------------------
% Initialize fit using SVD of ridge regression estimate
% ---------------------------------------------------
% do SVD on each relevant portion of w0
wridge = (XX + lam0*eye(nintot))\XY; % ridge regression solution 
wridgefilts = mat2cell(wridge,nin,nout); % convert into cell array for individual filters
wU = cell(ninpops,1);
wVt = cell(ninpops,1);
switch opts.Init
case 'rrr'
    for jj = 1:ninpops
        [~,~,vrrr] = svd(Yout'*Xin{jj}*wridgefilts{jj});  % perform SVD of relevant matrix
        vrrr = vrrr(:,1:rnks(jj));      % get column vectors
        % urrr= (X'*X)\(X'*Y*vrrr);  % get row vectors
        wU{jj}= wridge * vrrr;  % get row vectors
        s = vecnorm(wU{jj});
        wU{jj} = wU{jj} ./ s;
        wVt{jj} = (vrrr .* s)';
        s = s';
    end
    fprintf('rrr initialization\n');
case 'svd'
    for jj = 1:ninpops
        [u0,s{jj},v0] = svd(wridgefilts{jj},'econ','vector'); % SVD
        ii = 1:rnks(jj); % indices of singular vectors to keep
        wU{jj} = u0(:,ii) / sqrt(nin(jj));  % left singular vectors
        wVt{jj} = sqrt(nin(jj)) * diag(s{jj}(ii))*v0(:,ii)'; % right singular vectors (weighted by singular values)
    end
    fprintf('svd initialization\n');
case 'random' % randn initialization
    for jj = 1:ninpops
        [u0,s,v0] = svd(wridgefilts{jj},'econ','vector'); % SVD
        s = (norm(wridgefilts{jj},'fro') / rnks(jj)) * ones(rnks(jj),1);
        wU{jj} = randn(nin(1),rnks(jj));
        wVt{jj} = s(1) * randn(rnks(jj),nout); % right singular vectors (weighted by singular values)
    end
    fprintf('random initialization\n');
end

% ---------------------------------------------------
% Set up VI
% ---------------------------------------------------
% initiation
ELBO = -Inf;
dElBO = inf;
wchange = inf;

% hyperparameters
a_0 = .001;
b_0 = .001;
c_0 = .001;
d_0 = .001;


% calculate variational parameters that do not change
a_vi = a_0 + nout * nsamp / 2;
c_vi = c_0 + nout / 2;
b_vi = a_vi * var(Yout(:));
% if we're not fitting theta, override initiated s with true values
if ~if_fit_theta
    s = ard_params.s;
end
for ii = 1:numel(s)
    d_vi{ii} = c_vi .* (s{ii}(1:rnks(ii)).^2) / (nout* nin(ii));
end


% need better initialization
% variational parameters initialization
% size(cell2mat(wVt))
if if_fit_V
    V_vi_cell = wVt;
else
    V_vi_cell = ard_params.V_vi_cell;
end
V_vi = vec(cell2mat(wVt)); % vec Vt
if if_fit_U
    U_vi_cell = wU;
else
    U_vi_cell = ard_params.U_vi_cell;
end
U_vi = cell2mat(cellfun(@(x) vec(x), U_vi_cell, 'uni',0));

% % test vec to mat transformation
% tmp = vec2mats(U_vi,nin,rnks); % U 
% tmp = mat2cell(reshape(V_vi,sum(rnks),nout),rnks,nout); % Vt

SigmaU_vi = speye(nU);
SigmaV_vi = speye(nV);


% build S_i
% Ss = reshape(kron(eye(nout),Xin{1}),nsamp*nout,nin,nout); 
Pu = commutation(nintot,totrnks);
Pv = commutation(totrnks,nout);


% extra initialization to use
% Lmat_x = speye(nintot);
% Lmat_y = speye(nout);
% Lmat_r = speye(rnks(1));
Lmat_rx = speye(sum(prod_rnk_nin));
% plt_cmp_matrices({SPuPS,save.SPuPS_multi});ef

% components that will be used and updated in calculation
YY_vec = vec(Yout)'*vec(Yout);
to_save.ELBO = nan(opts.MaxIter,1);
iter = 1;
% ---------------------------------------------------
% Optimize:  Run alternating coordinate ascent on U and Vt
% ---------------------------------------------------
while (iter <= opts.MaxIter) && (abs(dElBO) > opts.TolFun)

    % saving old parameters
    % params_old = struct('U_vi',U_vi,'V_vi',V_vi,'SigmaU_vi',SigmaU_vi,'SigmaV_vi',SigmaV_vi,'a_vi',a_vi,'b_vi',b_vi,'c_vi',c_vi,'d_vi',d_vi);

    % update qV
    if if_fit_V
        sigma_mumu_U = SigmaU_vi + U_vi * U_vi';
        
        tmp = [];
        for ii = 1:numel(rnks)
            tmp = [tmp repmat(nin(ii),1,rnks(ii))];
        end
        blocks = mat2cell(sigma_mumu_U,tmp,tmp);
        % add those blocks to Au
        Au = sparse(nintot*totrnks, nintot*totrnks); 
        block_r = 0;
        for ii = 1:ninpops
            for mm = 1:rnks(ii)
                block_c = 0;
                for jj = 1:ninpops
                    for nn = 1:rnks(jj)
                        % block index
                        rs = block_r + sum(nin(1:ii-1)) + (1:nin(ii));
                        cs = block_c + sum(nin(1:jj-1)) + (1:nin(jj));
                        Au(rs,cs) = blocks{block_r/nintot + 1 , block_c/nintot + 1};
                        % fprintf('index Au: %d %d\n',block_r/nintot + 1 , block_c/nintot + 1);
                        % fprintf('mm, nn, ii, jj, r1, rend, c1, cend: %d %d %d %d %d %d %d %d\n',mm,nn,ii,jj,rs(1),rs(end),cs(1),cs(end));
                        block_c = block_c + nintot;
                    end
                end
                block_r = block_r + nintot;
            end
        end

        Au = mat2cell(Pu*Au*Pu', totrnks* ones(1,nintot),totrnks*ones(1,nintot));
        Au = arrayfun(@(xi,Aui) xi*Aui{1}, XX, Au, 'uni',0);
        SPuPS = kron(speye(nout),sum(cat(3,Au{:}),3));
        
        E_Theta = diag(repmat(c_vi ./ cat(1,d_vi{:}), nout,1));

        SigmaV_vi_inv = E_Theta + a_vi / b_vi * SPuPS;
        SigmaV_vi = inv(SigmaV_vi_inv);

        tmp = cellfun(@(u,x) u'*x',U_vi_cell,Xin','uni',0);
        V_vi = a_vi / b_vi * SigmaV_vi * vec(cell2mat(tmp)*Yout);
        V_vi_cell = mat2cell(reshape(V_vi,sum(rnks),nout),rnks,nout); % Vt
    end
    
    
    % update qU
    if if_fit_U
        sigma_mumu_V = SigmaV_vi + V_vi * V_vi';

        blocks = mat2cell(Pv*sigma_mumu_V*Pv',ones(1,totrnks)*nout,ones(1,totrnks)*nout);
        Av = cellfun(@(x) trace(x),blocks);

        % add those blocks to Av
        SPvPS = sparse(sum(nin.*rnks), sum(nin.*rnks)); 
        for ii = 1:ninpops % row index
            for jj = 1:ninpops % col index
                tmp = kron(Av((1:rnks(ii)) + sum(rnks(1:ii-1)), (1:rnks(jj)) + sum(rnks(1:jj-1))),...
                           XX((1:nin(ii)) + sum(nin(1:ii-1)),(1:nin(jj)) + sum(nin(1:jj-1))));
                SPvPS((1:prod_rnk_nin(ii)) + sum(prod_rnk_nin(1:ii-1)), (1:prod_rnk_nin(jj)) + sum(prod_rnk_nin(1:jj-1))) = tmp;
            end
        end
        
        % !!! Sigma_U_vi having negative values
        SigmaU_vi_inv = Lmat_rx + a_vi / b_vi * SPvPS;
        SigmaU_vi = inv(SigmaU_vi_inv);
        U_vi = a_vi / b_vi * SigmaU_vi * cell2mat(cellfun(@(x,v) vec(x'*Yout*v'), Xin', V_vi_cell, 'uni',0));
        U_vi_cell = vec2mats(U_vi,nin,rnks)';
    end


    % update other variational parameters
    tmp = cellfun(@(x,u,vt) x*u*vt, Xin', U_vi_cell, V_vi_cell, 'uni',0); % xuv';
    E_UV = YY_vec - 2*vec(Yout)'*vec(sum(cat(3,tmp{:}),3)) + trace(SPvPS*SigmaU_vi)+ U_vi'*SPvPS*U_vi;
    if if_fit_theta
        mumu_sig_V = reshape(diag(SigmaV_vi) + V_vi.^2,totrnks,[]);
        b_vi = b_0 + 1/2 * E_UV;
        d_vi = mat2cell(sum(mumu_sig_V,2)/2 + d_0,rnks)';
    end

    % calculate ELBO
    ELBO_old = ELBO;
    % ELBO = (nsamp*nout/2 + a_0 - a_vi) * psi(a_vi) + (totrnks*nout/2 + c_0 - c_vi) * psi(c_vi) ...
    %         - (nsamp*nout/2 + a_0) * log(b_vi) - (totrnks*nout/2 + c_0) * log(d_vi) ...
    %         - (E_UV/2 + b_0) * a_vi/b_vi - (mumu_sig_V/2 + d_0) * c_vi/d_vi ...
    %         - (U_vi'*U_vi + trace(SigmaU_vi))/2 + (my_logdet(SigmaU_vi) + my_logdet(SigmaV_vi))/2 ...
    %         + gammaln(a_vi) + a_vi + gammaln(c_vi) + c_vi;
    d_vi_all = cell2mat(d_vi');
    log_d_vi_all = log(d_vi_all);
    ELBO = nsamp*nout/2  * (psi(a_vi)-log(b_vi)) - a_vi / (2*b_vi) * E_UV ... % E[log p(Y|X,U,V,a,b)]
            - 1/ 2 * (U_vi'*U_vi + trace(SigmaU_vi)) ... % E[log p(U)]
            + nout / 2 * sum(psi(c_vi) - log_d_vi_all) - sum(repmat(c_vi./d_vi_all,nout,1) .* (diag(SigmaV_vi) + V_vi.^2)) / 2 ... % E[log p(V,c,d)]
            + (c_0 -1) * sum(psi(c_vi) - log_d_vi_all) - d_0 * sum(c_vi ./ d_vi_all) ... % E[log p(theta) p(c,d)]
            + (a_0 - 1) * (psi(a_vi) - log(b_vi)) - b_0 * a_vi / b_vi ... % E[log p(a,b)]
            + 0.5 * (my_logdet(SigmaU_vi) + my_logdet(SigmaV_vi)) ...
            - (log(b_vi) - gammaln(a_vi) - a_vi + (a_vi - 1) * psi(a_vi)) ...
            - sum(log_d_vi_all - gammaln(c_vi) - c_vi + (c_vi - 1) * psi(c_vi));
    % ELBO = nsamp*nout/2  * (-log(b_vi)) - a_vi / (2*b_vi) * E_UV + ...
    %         - 1/ 2 * (U_vi'*U_vi + trace(SigmaU_vi)) ...
    %         + nout / 2 * sum(psi(c_vi) - log(d_vi)) - sum(diag(c_vi ./ d_vi) * mumu_sig_V,"all") / 2 ...
    %         + (c_0 -1) * sum(psi(c_vi) - log(d_vi)) - d_0 * sum(c_vi ./ d_vi) ...
    %         + (a_0 - 1) * ( - log(b_vi)) - b_0 * a_vi / b_vi ...
    %         + 0.5 * (my_logdet(SigmaU_vi) + my_logdet(SigmaV_vi)) ...
    %         - (log(b_vi) - gammaln(a_vi) - a_vi) ...
    %         - sum(log(d_vi) - gammaln(c_vi) - c_vi + (c_vi - 1) * psi(c_vi));
    % ELBO2 = psi(a_vi) * (nsamp*nout/2 + a_0 - a_vi)
    % ELBO - ELBO2

    
    dElBO = ELBO - ELBO_old;
    to_save.ELBO(iter) = ELBO;
    % saving new parameters
    % params_new = struct('U_vi',U_vi,'V_vi',V_vi,'SigmaU_vi',SigmaU_vi,'SigmaV_vi',SigmaV_vi,'a_vi',a_vi,'b_vi',b_vi,'c_vi',c_vi,'d_vi',d_vi);
    
    % bp;
    if mod(iter,100)==0
        fprintf('iter %d  ELBO: %f, dELBO: %f\n',iter,ELBO,dElBO);
    end
    iter = iter + 1;
    % tmp_plt;
end


% check parameters
% plt_cmp_matpairs([U_vi_cell,ard_params.U_vi_cell]);export_fig('tmp.png','-m3'); % not that accurate
% plt_cmp_matpairs([V_vi_cell,ard_params.V_vi_cell]);export_fig('tmp.png','-m3'); % not that accurate
% plt_cmp_matpairs([cellfun(@(u,vt) u*vt, U_vi_cell, V_vi_cell, 'uni',0) cellfun(@(u,vt) u*vt, ard_params.U_vi_cell, ard_params.V_vi_cell, 'uni',0)]);export_fig('tmp.png','-m3'); 
% cellfun(@(u,vt) u*vt, ard_params.U_vi_cell, ard_params.V_vi_cell, 'uni',0);
% ard_params.U_vi_cell


to_save.SigmaU_vi = SigmaU_vi;
to_save.SigmaV_vi = SigmaV_vi;
to_save.a_vi = a_vi;
to_save.b_vi = b_vi;
to_save.c_vi = c_vi;
to_save.d_vi = d_vi;


% reverse vectorization
% U_vi = reshape(U_vi,nin(1),rnks(1));
% V_vi = reshape(V_vi,rnks(1),nout);
% B = U_vi * V_vi;

% save tmp;error;
% to_save.niter = iter-1;
% to_save.lambda = alpha * nsevar;

% if (fchange > opts.TolFun)
%     warning('bilinearMultifiltRRR_coordAscent: did not converge to desired tolerance');
% end

end
