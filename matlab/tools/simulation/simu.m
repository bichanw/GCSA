classdef simu < handle
% create simulation
	methods (Static)

        function cmp_ARD()
            % compare different ARD methods, which gives you better resutls
            [X_fit,Y_hist_fit,Y_fit,tr_ind,flts,ops.sim_params] = simu.AR1(struct('rnk',[5 4],'magnitude',[1 1],'signse',0.5,'nt',5,'ntr',100));


            Xin  = {X_fit,Y_hist_fit{:}};
            Yout = Y_fit;
            rnks = ops.sim_params.rnk;
            lam0 = 10;
            opts = struct('MaxIter',1e3,'TolFun',1e-6,'Display','iter');
            
            % fixed lambda
            [wU,wVt,wwfilts,to_save] = bilinearMultifiltRRRRidgeARD_coordAscent(Xin,Yout,rnks,lam0,opts);

            % compare with old code
                % % fixed lambda, correct code
                % [wU(:,end+1),wVt(:,end+1),wwfilts(:,end+1),to_sav2] = bilinearMultifiltRRR_coordAscent(Xin,Yout,rnks,lam0,opts);
                % % compare results
                % plt_cmp_matrices({wwfilts{1,1},wwfilts{1,2},flts(1).wfilt},{'ortho','old','true'});ef;


            % test if U/V are both semi-orthogonal
            ax = np(1,2);
            myimg(ax(1),wU{1}'*wU{1});
            myimg(ax(2),wVt{1}*wVt{1}');
            title(ax(1),'U^T U'); title(ax(2),'V^T V');
            arrayfun(@(h) colorbar(h),ax);
            ef;


            % ARD finding the best lambda, then use it for parameter fitting
            [kridge,hprs_hat] = autoRidgeEfficient_fixedpoint(cell2mat(Xin),Yout,1e3,struct('maxiter',1e3,'tol',1e-6));
            hprs_hat.alpha * hprs_hat.nsevar
            [wtfit,wxfit,wwfilts2,to_save] = bilinearMultifiltRRR_coordAscent(Xin,Yout,rnks,hprs_hat.alpha * hprs_hat.nsevar,struct('Display','none','MaxIter',1e3));

            plt_cmp_matrices({wwfilts{1},wwfilts2{1},flts(1).wfilt},{'joint','ARD first','true'});ef;
            % plt_cmp_matrices({wU{1},wtfit{1},flts(1).wx},{'joint','ARD first','true'});ef;
            % plt_cmp_matrices({wVt{1},wxfit{1},flts(1).wy},{'joint','ARD first','true'});ef;
        end

        function rank2_vs_copies()
            % still using AR1 model
            [X_fit,Y_hist_fit,Y_fit,tr_ind,flts,ops.sim_params] = simu.AR1(struct('rnk',[1 5],'magnitude',[1 1],'signse',0.5,'nt',5,'ntr',500));
        
            % parameters
            Xin  = {X_fit,Y_hist_fit{:}};
            Yout = Y_fit;
            rnks = ops.sim_params.rnk;
            lam0 = 10;
            opts = struct('MaxIter',1e3,'TolFun',1e-6,'Display','iter');



            % regularization strength
            [kridge,hprs_hat] = autoRidgeEfficient_fixedpoint(cell2mat({X_fit,Y_hist_fit{1}}),Y_fit,10,struct('maxiter',1e3,'tol',1e-6));
            lambda = hprs_hat.lambda;
            
            % rank 2 model
            [B,ops_] = regressors.bilinear.train({X_fit,X_fit,Y_hist_fit{1}},Y_fit,[2 0 9],struct('lambda',lambda));
            % 2 copies model
            [B(end+1,:),ops_(end+1)] = regressors.bilinear.train({X_fit,X_fit,Y_hist_fit{1}},Y_fit,[1 1 9],struct('lambda',lambda));

            % compare the filters
                ax = plt_cmp_matrices({B{1,1},B{2,1}+B{2,2}},{'rank 2','2 copies of rank 1'});
                title(ax(3),'\Delta','Interpreter','tex');
                align_ax(ax([1 4]),false,false,true);
                set(gcf,'Position',[0 0 300 225]);
                ef;

            % U / V filters for 2 copies
            % without doing svd
                ax = np(2,2);
                plt_cmp_vectors(ops_(1).wtfit{1},{'flt 1','flt 2'},'lines',ax(1));
                plt_cmp_vectors(ops_(1).wxfit{1}',{'flt 1','flt 2'},'lines',ax(3));
                plt_cmp_vectors([ops_(2).wtfit{1}, ops_(2).wtfit{2}],{'flt 1','flt 2'},'lines',ax(2));
                plt_cmp_vectors([ops_(2).wxfit{1}; ops_(2).wxfit{2}]',{'flt 1','flt 2'},'lines',ax(4));
                t = text_legend(ax(1),{'flt 1','flt 2'},[],1);
                title(ax(1),'rank 2'); title(ax(2),'2 copies');
                ylabel(ax(1),'U'); ylabel(ax(3),'V');
                ef;

            % higher rank?
            [U,S,V] = svd(B{1,1},"vector","econ");
        end


        

        function main()
            % run simulation multiple times and see prediction profile

            for ii = 1:10
                % create simulation
                [X_fit,Y_hist_fit,Y_fit,tr_ind,flts,ops.sim_params] = simu.AR1(struct('rnk',[5 4],'magnitude',[1 1],'signse',0.1));
                [B,ops_cv] = regressors.bilinear.train({X_fit,Y_hist_fit{1}},Y_fit,[5 4],struct('lambda',0));
                plt_cmp_matrices({B{1}, flts(1).wfilt});ef;

                % [X_fit,Y_hist_fit,Y_fit,tr_ind,flts,ops.sim_params] = simu.recreate_grid(struct('signse',3,'Y2Xnse',10,'sig_U',0));
                % [X_fit,Y_hist_fit,Y_fit,tr_ind,flts,ops.sim_params] = simu.common_input(struct('sig_U',2,'magnitude',[ones(2,1)*1 ones(2,1)*0.8]));
                % [X_fit,Y_hist_fit,Y_fit,tr_ind,flts,ops.sim_params,raw_data] = simu.entangled(struct('ntr',100,'nt',10,'signse',5,'magnitude',[1*ones(2,1) 2*ones(2,1)],'rnk',[4 4; 4 4],'seed',3));
                % ax = np(2,2); for ii = 1:4 imagesc(ax(ii),flts(ii).wfilt);colorbar(ax(ii)); end % check filter magnitude

                % fit using biliRRR script
                ops.spk_count = struct('tr_ind',tr_ind);
                ops.Ranks_x = 0:6;
                ops.Ranks_yhist = 0:6;
                ops.n_fold = 5;
                % [Err,ops_cv] = plt.cv_grid(X_fit,Y_hist_fit,Y_fit,ops,'',true);
                [X_fit,Y_fit,Y_hist_fit,tr_ind] = simu.organize(raw_data.X(:,2:end,:),raw_data.Y);
                [Err,ops_cv] = plt.cv_grid(X_fit,Y_hist_fit,Y_fit,ops,sprintf('sim_entangled_%d',ii),true);
                [X_fit,Y_fit,Y_hist_fit,tr_ind] = simu.organize(raw_data.X(:,1:end-1,:),raw_data.Y);
                [Err,ops_cv] = plt.cv_grid(X_fit,Y_hist_fit,Y_fit,ops,sprintf('sim_entangledt-1_%d',ii),true);
                [X_fit,Y_fit,Y_hist_fit,tr_ind] = simu.organize(raw_data.Y(:,2:end,:),raw_data.X);
                [Err,ops_cv] = plt.cv_grid(X_fit,Y_hist_fit,Y_fit,ops,sprintf('sim_entangledreverse_%d',ii),true);

            end
   
        end

        function from_data()
            % extract X, Y0, B from data and simulate, recover
            ops = struct('n_fold',10,'lambda',1e3,'Ranks_x',0:9,'Ranks_yhist',[0 9],'rnk',flip(combinations([0 9],0:9),2)); % saving ops
            [spk,ops.spk_count,X] = count_spk.multi_bins(data,data.events.feedback_times,{'MO'},struct('bin_width',0.1,'n_bins',1/0.1,'n_hist',1,'shuffle_trials',false,'subtract_PSTH',false));
            [spk(end+1),~,Y] = count_spk.multi_bins(data,data.events.feedback_times,{'DLS'},struct('bin_width',0.1,'n_bins',1/0.1,'n_hist',1,'shuffle_trials',false,'subtract_PSTH',false));
        
            % fit Bx, Bt-1
            [B,ops.train] = regressors.bilinear.train({squeeze(spk{1}(1,:,:)),squeeze(spk{2}(2,:,:))},squeeze(spk{2}(1,:,:)),[1 9]);
            
            % remove cell in X, Y with 0 variance
            % X = X(:,:,~ops.train.exclude_id{1});
            % Y = Y(:,:,~ops.train.exclude_id{2});
            % reshape(ops.spk_count.tr_ind,[],numel(unique(ops.spk_count.tr_ind)))
            % tmp = reshape(squeeze(spk{2}(2,:,:)),[],numel(unique(ops.spk_count.tr_ind)),size(spk{2},3));
            Y0 = squeeze(Y(1,:,:));

            % simulate data from B
            Y_sim = nan(size(Y));
            Y_sim(1,:,:) = Y0;
            signse = 0.1;
            for ii = 2:size(Y,1)
                Y_sim(ii,:,:) = squeeze(Y_sim(ii-1,:,~ops.train.exclude_id{2})) * B{2} + squeeze(X(ii,:,~ops.train.exclude_id{1})) * B{1} + signse*randn(size(Y,2),size(Y,3));
            end

            % can we recover B from the simulated data?
            spk_simed = count_spk.spk_mat2cell(Y_sim,1);
            [B_sim,ops.train_sim] = regressors.bilinear.train({squeeze(spk{1}(1,:,:)),squeeze(spk_simed(2,:,:))},squeeze(spk_simed(1,:,:)),[1 9]);
            % compare B and B_sim
            plt_cmp_matrices({B{1},B_sim{1}(~ops.train.exclude_id{1},:)},{'fitted','sim'}); ef;%export_fig('B_x with sim.pdf');
            plt_cmp_matrices({B{2},B_sim{2}(~ops.train.exclude_id{2},:)},{'fitted','sim'}); ef;%export_fig('B_yhist with sim.pdf');

            % plot grid
                ops = rmfield(ops,'Xmean');
                [err,ops,to_save] = myCV.by_trial({squeeze(spk{1}(1,:,:)),squeeze(spk_simed(2,:,:))},squeeze(spk_simed(1,:,:)),ops.spk_count.tr_ind,ops);
                % plt_CSA.grid(err,ops.rnk);
                plt_CSA.predicatibility(reshape(err,10,2,[]),0:9,np);
                clear tmp; tmp(1,:,:) = err;
                ndim = loss.cal_rnk(tmp(:,ops.rnk(:,2)==0,:),0:9);
                ndim(2) = loss.cal_rnk(tmp(:,ops.rnk(:,1)==1,:),0:9)

                ef;

            % are the activities similar?
            ax = plt_cmp_matrices({squeeze(Y(2,:,:)),squeeze(Y_sim(2,:,:))},{'true','sim'});
            align_ax(ax([1 4]),false,false,true);
            Y_true = squeeze(spk{2}(end,:,:));
            Y_pred = squeeze(spk_simed(end,:,:));
            tr_ind = ops.spk_count.tr_ind;
            ax = plt.cmp_y_pred(squeeze(spk{2}(end,:,:)),squeeze(spk_simed(end,:,:)),ops.spk_count.tr_ind,5);
            % what about adding within sample prediction
            Y_pred = regressors.bilinear.predict({squeeze(spk{1}(end,:,:)),squeeze(spk{2}(end-1,:,:))},B,ops.train);
            ax = plt.cmp_y_pred(squeeze(spk{2}(end,:,:)),{squeeze(spk_simed(end,:,:)),Y_pred},ops.spk_count.tr_ind,5);
        end

        function [X,Y,B,ops] = linear(ops)
            % simplest case: Y = XB
            if nargin == 0
                ops = struct;
            end
            [T,ops] = getOr(ops,'T',200);
            [nx,ops] = getOr(ops,'nx',10);
            [ny,ops] = getOr(ops,'ny',5);
            [signse,ops] = getOr(ops,'signse',0.1);
            [magnitude,ops] = getOr(ops,'magnitude',1);
            
            X = randn(T,nx);
            B = randn(nx,ny) * magnitude;
            Y = X * B + signse*randn(T,ny);
        end

        function [X,Y,U,V,ops] = RRR(ops)
            % simplest case: Y = XB
            if nargin == 0
                ops = struct;
            end
            [T,ops] = getOr(ops,'T',100);
            [nx,ops] = getOr(ops,'nx',10);
            [ny,ops] = getOr(ops,'ny',5);
            [rnk,ops] = getOr(ops,'rnk',1);
            [signse,ops] = getOr(ops,'signse',0.1);
            [magnitude,ops] = getOr(ops,'magnitude',1);
            [thetas,ops] = getOr(ops,'thetas',[]);
            [Sigma,ops] = getOr(ops,'Sigma',[]);
            % input U, V if not wish to generate a new pair
            [U,ops] = getOr(ops,'U',[]);
            [V,ops] = getOr(ops,'V',[]);
            

            % generate input
            X = randn(T,nx);
            X = X - mean(X,1); % add an extra step to zero data?

            % generate communication matrix
            if isempty(U)
                if isempty(thetas)
                    U = randn(nx,rnk);
                    V = randn(ny,rnk);
                    [U,S,V] = svd(U*V','econ');
                    U = U(:,1:rnk);
                    V = V(:,1:rnk);

                    % scale V by theta input (for ARD we assume U as orthonormal)
                    U = U*S(1:rnk,1:rnk);
                    V = V * magnitude;
                else % generate in a way that fits ARD assumption
                    U = randn(nx,rnk);
                    V = randn(ny,rnk) ./ sqrt(thetas);
                end
            end
            W = U*V';

            % generate Y
            if isempty(Sigma)
                E = signse*randn(T,ny);
                % E = E - mean(E,1);
                Y = X * W + E;
            else
                Y = X * W + mvnrnd(zeros(1,ny),Sigma,T);
            end
            % tmp  =mvnrnd(zeros(1,ny),Sigma,T);
        end
        
        function [X_fit,Y_hist_fit,Y_fit,tr_ind,flts,params] = independent_U(params)
            % try to recreate grid pattern from data

            % simulation 1:
            % common input source U to X & Y
            % recurrent connections within Y

            % X = U * Bx; Bx is full rank
            % Y = U * By + Yt-1 * Byhist; both By and Byhist are low rank

            % parameter intitiation
            if nargin < 1
                params = struct;
            end
            [ntr,params] = getOr(params,'ntr',100);
            [nt ,params] = getOr(params,'nt',5);
            [nx ,params] = getOr(params,'nx',8);
            [ny ,params] = getOr(params,'ny',12);
            [nu ,params] = getOr(params,'nu',8);
            [rnk,params] = getOr(params,'rnk',[4,3]);  % rank of each coefficient matrix; vector of 2 for AR1
            [signse,params] = getOr(params,'signse',0.3);
            [U2Xnse,params] = getOr(params,'U2Xnse',0.1);
            nhist = 1; % AR1 for Y
            

            % generate coefficient
            flts = gen_flts([nu, nx; nu,ny; ny,ny],[min([nu nx]), rnk]);

            % -------------------------------------
            % preventing blowing up for history effect
            B = concate_filters(flts(2:end));
            scale = 1; 
            alpha = max(abs(eig(B)));
            while alpha > 1
                fprintf('alpha = %.2f\n',alpha);
                % scale coefficients and noise mangitude
                for jj = 1:numel(flts)
                    flts(jj).wfilt = flts(jj).wfilt / (alpha+ 0.1);
                end
                scale = scale / (alpha + 0.1);
                B = concate_filters(flts(2:end));
                alpha = max(abs(eig(B)));   
            end
            signse = signse * scale;
            U2Xnse = U2Xnse * scale;
            params.signse = signse;
            params.U2Xnse = U2Xnse;


            % -------------------------------------
            % generate intial data
            U = randn(ntr,nt,nu);
            Y = randn(ntr,nhist,ny);
            X = nan(ntr,nt,nx);

            % generate full X and Y
            for jj = 1:nt

                % random noise
                Y(:,nhist+jj,:) = signse*randn(ntr,ny);
                X(:,jj,:) = squeeze(U(:,jj,:)) * flts(1).wfilt + U2Xnse*randn(ntr,nx);
                
                % input  
                Y(:,nhist+jj,:) = squeeze(Y(:,nhist+jj,:)) + squeeze(U(:,jj,:)) * flts(2).wfilt;

                % own histogry effect
                for kk = 1:nhist
                    Y(:,nhist+jj,:) = squeeze(Y(:,nhist+jj,:)) + squeeze(Y(:,nhist+jj-kk,:)) * flts(3).wfilt;
                end
            end


            % organize data for fitting
            [X_fit,Y_fit,Y_hist_fit,tr_ind] = simu.organize(X,Y);


        end


        function [X_fit,Y_hist_fit,Y_fit,tr_ind,flts,params] = AR1_jp(params,flts)

            % put the code here to figure out why results were different from jonathan's code
            % compare filters: 
            nin = 8; % number of neurons in intput population
            nout = 10; % number of neurons in output population
            ntr = 50;
            rnk1 = 1;  % rank of input weights 
            rnk2 = 5;  % rank of auto-regressive weights


            % generate filters with my code
            flts = gen_flts([nin nout; nout nout],[rnk1 rnk2],[1 5]);
            [flts,scale] = scale_flts(flts);

            ax = np(2);
            arrayfun(@(ii) myimg(ax(ii),flts(ii).wfilt),1:2);
            arrayfun(@(ii) colorbar(ax(ii)),1:2);
            ef;

            % generate filters with jonathan's code
            flts = gen_flts_jp(nin,nout,ntr,rnk1,rnk2);

        end


        function [X_fit,Y_hist_fit,Y_fit,tr_ind,flts,params] = AR1(params,flts)
            % simulate data baesd on AR1 assumption

            % parameter intitiation
            if nargin < 1
                params = struct;
            end
            [ntr,params] = getOr(params,'ntr',200);
            [nt ,params] = getOr(params,'nt',1);
            [nx ,params] = getOr(params,'nx',12);
            [ny ,params] = getOr(params,'ny',15);
            [self_couple, params] = getOr(params,'self_couple',0);
            [rnk,params] = getOr(params,'rnk',[2,1]);  % rank of each coefficient matrix; vector of 2 for AR1
            [signse,params] = getOr(params,'signse',0.8);
            nhist = numel(rnk)-1; % how many history points to use

            [magnitude,params] = getOr(params,'magnitude',[5,1*ones(1,nhist)]); % introduce different magnitudes for each filter
            params.seed = rng('Shuffle');

            % -------------------------------------
            if nargin < 2 || isempty(flts)

                % generate filter 
                flts = gen_flts([nx ny; ny ny],rnk,magnitude); fprintf('use my code gen\n');
                % flts = gen_flts_jp(nx,ny,ntr,rnk(1),rnk(2)); fprintf('use jp gen\n');

                % add self-coupling to the last filter
                flts(end).wfilt = flts(end).wfilt + eye(ny) * params.self_couple;
            end

            % -------------------------------------
            % preventing blowing up for history effect
            B = concate_filters(flts);
            scale = 1; 
            alpha = max(abs(eig(B)));
            while alpha > 1
                fprintf('alpha = %.2f\n',alpha);
                % scale coefficients and noise mangitude
                for jj = 1:numel(rnk)
                    flts(jj).wfilt = flts(jj).wfilt / (alpha / 0.98);
                end
                scale = scale / (alpha  / 0.98);
                B = concate_filters(flts);
                alpha = max(abs(eig(B)));   
            end
            % [flts,scale] = scale_flts(flts);
            signse = signse * scale;
            params.signse = signse;

            % -------------------------------------
            % generate random X and Y
            X = randn(ntr,nt,nx) * scale;
            Y = randn(ntr,nhist,ny) * scale * signse; % starting history of Y

            % smooth X on the second dimension
            for ii = 1:size(X,3)
                X(:,:,ii) = gsmooth(squeeze(X(:,:,ii)),3);
            end
            % function g = gsmooth(y, sig)
            % gsmooth - smooth vector or matrix by filtering with a Gaussian 
            %
            % g = gsmooth(x, sig);
            %
            % Inputs:
            %     y [MxN] - matrix or vector (if matrix, operates along columns only)  
            %   sig [1x1] - stdev of smoothing Gaussian (in samples)

            % generate full Y
            if size(Y,1) == 1

                X = squeeze(X);
                Y = squeeze(Y)';
                % allocate Y space
                Y = cat(1,Y,nan(nt,ny));

                for jj = 1:nt
                    
                    % random noise and input
                    Y(nhist+jj,:) = signse*randn(ntr,ny) + X(jj,:) * flts(1).wfilt;
                    
                    % own histogry effect
                    for kk = 1:nhist
                        Y(nhist+jj,:) = Y(nhist+jj,:) + Y(nhist+jj-kk,:) * flts(kk+1).wfilt;
                    end


                end
            else
                for jj = 1:nt
                        
                    % random noise
                    Y(:,nhist+jj,:) = signse*randn(ntr,ny);
                    
                    if size(Y,1) == 1
                        Y(:,nhist+jj,:) = Y(:,nhist+jj,:) + X(:,jj,:) * flts(1).wfilt;
                        Y(:,nhist+jj,:) = squeeze(Y(:,nhist+jj,:))' + squeeze(X(:,jj,:))' * flts(1).wfilt;
                        % own histogry effect
                        for kk = 1:nhist
                            Y(:,nhist+jj,:) = squeeze(Y(:,nhist+jj,:))' + squeeze(Y(:,nhist+jj-kk,:))' * flts(kk+1).wfilt;
                        end
                    else
                        % input  
                        Y(:,nhist+jj,:) = squeeze(Y(:,nhist+jj,:)) + squeeze(X(:,jj,:)) * flts(1).wfilt;
                        % own histogry effect
                        for kk = 1:nhist
                            Y(:,nhist+jj,:) = squeeze(Y(:,nhist+jj,:)) + squeeze(Y(:,nhist+jj-kk,:)) * flts(kk+1).wfilt;
                        end
                    end


                end
            end

            % -------------------------------------
            % organize data for fitting
            if ntr == 1
                X_fit = X;
                Y_fit = Y(2:end,:);
                Y_hist_fit = {Y(1:end-1,:)};
                tr_ind = [];
            else
                [X_fit,Y_fit,Y_hist_fit,tr_ind] = simu.organize(X,Y);
            end

        end

            
        % organize data for fitting process
        function [X_fit,Y_fit,Y_hist_fit,tr_ind] = organize(X,Y)
            % find parameters from input
            % resulting X, Y should be samples * cells
            [ntr,nt,nx] = size(X);
            [~,tmp,ny]  = size(Y); nhist = tmp - nt;

            X_fit = [];
            Y_fit = [];
            for jj = 1:nt
                X_fit = cat(1,X_fit,squeeze(X(:,jj,:)));
                Y_fit = cat(1,Y_fit,squeeze(Y(:,jj+nhist,:)));
            end
            Y_hist_fit = cell(1);
            for ii = 1:nhist
                Y_hist_fit{ii} = [];
                for jj = 1:nt
                    Y_hist_fit{ii} = cat(1,Y_hist_fit{ii},squeeze(Y(:,jj-ii+nhist,:)));
                end
            end
            tr_ind = ones(nt,1) * (1:ntr); % trial index
            tr_ind = tr_ind(:);
        end

        
    end

end