% add a path definition to automatically load?


%% technical paper
% fig 4: ARD - one input
    % figure /- compare errors and thetas?
        % params.thetas
        % ax = np;plot(to_save.c_vi ./ to_save.d_vi{1} );ef; % 3 5 50
        % ax = np;plot( to_save.d_vi{1} ./ to_save.c_vi);ef; % 3 5 50
        % load('results/ard_figure.mat');

        close all;
        ax = arrayfun(@(ii) subplot(1,2,ii,'Nextplot','add'), 1:2);
        h = mybar(Errs*100,ax(1));
        set(h.lines,'Color',[0 0 0 0.2]);
        set(ax(1),'XTick',[1 2],'XTickLabel',{'w/ ARD','RRR'},...
               'XLim',[0.5 2.5],'YLim',[0 0.2],'YTick',0:0.1:0.2);
        ylabel(ax(1),'mean-squared error (%)');
        % ax = np;
        plot_multiple_lines(c_vi ./ d_vis,ax(2),'x',1:5);
        % plot_multiple_lines( d_vis / c_vi,ax(2),'x',1:5);
        set(ax(2),'XLim',[0.5 5.5],'XTick',1:5);
        xlabel(ax(2),'i-th communication axis');
        ylabel('c_i / d_i');
        set(ax,'FontSize',8);
        % change font size of labels to 9
        set(gcf,'Position',[0 0 4 1.7]*72);
        exportgraphics(gcf,'tmp.png','Resolution',300);
        
    return

    % run and fit
    nsims = 100;
    Errs  = nan(nsims,2);
    rnks = 5; d_vis = nan(nsims,rnks);
    for isim = 1:nsims
        switch 1
        case 1
            % rrr simu
            [X_fit,Yout,U,V,params] = simu.RRR(struct('thetas',[5 .5 .5],'signse',1,'magnitude',1,'nx',10,'ny',8,'rnk',3,'T',50));
            flts = struct('wfilt',U*V');
            Xin = {X_fit};
        case 2
            % ar1 simu
            [X_fit,Y_hist_fit,Yout,tr_ind,flts,params] = simu.AR1(struct('rnk',[2 3],'signse',.5,'ntr',30,'nt',3));
            Xin = {X_fit,Y_hist_fit{1}};
        end
        lam0 = 0;
        opts = struct('MaxIter',1e5,'TolFun',1e-6,'Display','iter','Init','svd');
        ard_params = params2ard_true(flts,params.rnk); 
        [U_vi_cell,V_vi_cell,to_save] = rankr_VI_multi(Xin,Yout,rnks,lam0,opts,ard_params);
            
        % calculate error for figures
        [B,ops_rrr] = bilinear.train(Xin,Yout,rnks,struct());
        % plt_cmp_matrices({flts(1).wfilt,U_vi_cell{1}*V_vi_cell{1},B{1}},{'true','ARD','rrr'});ef;
        err_fun = @(a,b) mean((a-b).^2,'all') / mean(b.^2,'all');
        fprintf('ard err %.3f, rrr err %.3f\n',err_fun(flts(1).wfilt,U_vi_cell{1}*V_vi_cell{1}),err_fun(flts(1).wfilt,B{1}));
        Errs(isim,1) = err_fun(U_vi_cell{1}*V_vi_cell{1},flts(1).wfilt);
        Errs(isim,2) = err_fun(B{1},flts(1).wfilt);
        d_vis(isim,:) = to_save.d_vi{1};
    end
    c_vi = to_save.c_vi;
    save('results/ard_figure2.mat','Errs','d_vis','c_vi');
    
    return

    % thinking one input works better, but ard gives similar range of results than rrr?
    % try to change communication channel strength?


% supp fig: with self-recurrence
    [f,Hemis,ses_w_file] = load_data.results('MO','DLS','cv_grid_wself_0and9');
    err = permute(cat(3,f(:).err),[3 1 2]);
    
    plt_CSA.rnk_scatter(err,f(1).ops.rnk);
    
    % compare history rank w/ or w/o 
        % M2: rank 9
        switch 2
        case 1
            [f,Hemis,ses_w_file] = load_data.results('MO','DLS','cv_grid_wself_vary_hist');
            [f2,Hemis,ses_w_file] = load_data.results('MO','DLS','cv_grid_vary_hist');
            err  = permute(cat(3,f(:).err),[3 1 2]);
            ndim = loss.cal_rnk(err(:,f(1).ops.rnk(:,1)==9,:),f(1).ops.rnk(f(1).ops.rnk(:,1)==9,2));
            err  = permute(cat(4,f2(:).Err),[4 1 2 3]);
            ndim2 = loss.cal_rnk(squeeze(err(:,:,2,:)),0:9);
        case 2 % same for input
            [f,Hemis,ses_w_file] = load_data.results('MO','DLS','cv_grid_wself_0and9');
            [f2,Hemis,ses_w_file] = load_data.results('MO','DLS','cv_grid_0and9');
            err  = permute(cat(3,f(:).err),[3 1 2]);
            ndim = loss.cal_rnk(err(:,f(1).ops.rnk(:,2)==9,:),f(1).ops.rnk(f(1).ops.rnk(:,2)==9,1));
            err  = permute(cat(4,f2(:).Err),[4 1 2 3]);
            ndim2 = loss.cal_rnk(squeeze(err(:,:,2,:)),0:9);
        end

        ax = np;
        scatter(ax,ndim2+randn(size(ndim2))*0.1,ndim+randn(size(ndim))*0.1);
        xlabel('w/o self recurrence');
        ylabel('w/ self recurrence');
        equalize_plot(ax,true);
        set(gcf,'Position',[0 0 150 120]);
        ef;

    % only hist, no recur
        % rank determination
        [f,Hemis,ses_w_file] = load_data.results('MO','DLS','only_hist_recur_only_hist_recur');
        err = permute(cat(3,f(:).err),[3 1 2]);
        plt_CSA.rnk_scatter(err,flip(combinations([0 15],0:15),2));ef;
    % test w/o recurrence entirely
        [f,Hemis,ses_w_file] = load_data.results('MO','DLS','only_hist_recur_only_hist_recur');
        [f2,Hemis,ses_w_file] = load_data.results('MO','DLS','cv_grid_wself_vary_hist');
        plt_CSA.predicatibility()
        tmp = flip(combinations([0 9],0:9),2);
        
        isession = 1;
        ax = np;
        % low_rank history, w/ w/o recurrence
        plt_CSA.predicatibility(reshape(f(isession).err,16,2,[]), 0:15, ax(1));ef;
        % low-rank history + input
        
        h = plot_multiple_lines(1-f2(isession).err(f2(isession).ops.rnk(:,1)==9,:)',ax(1));ef;
        % h = plot_multiple_lines(data,ax,varargin)
        % plot the average lines while keeping individual traces
        % rows are individiual samples
        % data: 		Nsample * Ntime 
        % parameters:	base_color
        
    % sanity check, cell count in DLS
        sessions_oi = plt.data_exam.examine_cellcount_samehemi({'MO','DLS'},16);
        regionnames = readcell(sprintf('data/cellcount%d.csv',0));
        regionnames = regionnames(2:end,1);
        ind = find(strcmpi(regionnames, 'DLS'));
        cellcount = readmatrix('data/cellcount0.csv');
        cellcount(:,:,end+1) = readmatrix('data/cellcount1.csv');
        ax = np; histogram(max(cellcount(ind,sessions_oi+1,:),[],3),'BinEdges',0:25:150);ef;



% fig 3: rank hist vs. input
    [f,Hemis,ses_w_file] = load_data.results('MO','DLS','cv_grid_0and9');
    [f(:,end+1),Hemis,ses_w_file] = load_data.results('MO','DLS','cv_grid_vary_hist');

    ndim2 = [];
    for ii = 1:2
        Err = permute(cat(4,f(:,ii).Err),[4 1 2 3]);
        ndim = loss.cal_rnk(squeeze(Err(:,:,1,:)),f(1).ops.Ranks_x);
        ndim(end+1,:) = loss.cal_rnk(squeeze(Err(:,:,end,:)),f(1).ops.Ranks_x);
        ndim2(ii,:,:) = ndim';
    end

    % ndim2 = reshape(ndim,15,2,[]);
    % plot shuffled scatter
        ndim = squeeze(ndim2(2,:,:));
        Colors = cbrewer2('Set1',2);
        ax = np;
        % shuffle off-diagonal points
        ind = ndim(:,1)==ndim(:,2);
        ind = false(size(ndim,1),1);
        tmp = ndim;
        tmp(~ind,:) = tmp(~ind,:) + randn(sum(~ind),2)*0.1;
        s = scatter(tmp(:,1),tmp(:,2));
        lim = [0.5 max(ndim(:))+0.5];
        plot(ax,lim,lim,'k--');
        set(ax,'XLim',lim,'YLim',lim);
        set(s,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',0.3,'SizeData',50); 
        xlabel(ax,'identified dim CSA','Color',Colors(1,:));
        ylabel(ax,{'identified dim','Recurrent CSA'},'Color',Colors(2,:));
        set(gcf,'Position',[0 0 280 240]);
        set(ax,'FontSize',18);
        ef;
    


% fig 3: rank adding NAc
    [f,Hemis,ses_w_file,Sessions] = load_data.results('MO','DLS','cv_grid','NAc');
    Err = cat(4,f(:).Err);
    Err = permute(Err,[4 1 2 3]);

    ax = np(numel(f));
    arrayfun(@(ii) plt_CSA.predicatibility(squeeze(Err(ii,:,:,:)), f(ii).ops.Ranks_x, ax(ii)), 1:numel(f));

% fig 3: x-y dim plot
    switch 2
    case 1
        [f,SessionNames] = load_data.semedo_results('cv_grid_sepstim_semedo');
        Err = cat(1,f(:).Err);
    case 2
        [f,Hemis,ses_w_file] = load_data.results('MO','DLS','cv_grid_vary_hist');
        % [f,Hemis,ses_w_file] = load_data.results('MO','DLS','cv_grid_0and9');
        Err = permute(cat(4,f(:).Err),[4 1 2 3]);
    end

    % average across dataset
    ops = f(1).ops;
    Colors = cbrewer2('Set1',4);
    
    % err: ntrials x nranks x nfold
    ndim = loss.cal_rnk(squeeze(Err(:,:,1,:)),ops.Ranks_x);
    ndim(end+1,:) = loss.cal_rnk(squeeze(Err(:,:,end,:)),ops.Ranks_x);
    ndim = ndim';

    % plot shuffled scatter
        Colors = cbrewer2('Set1',2);
        ax = np;
        % shuffle off-diagonal points
        ind = ndim(:,1)==ndim(:,2);
        ind = false(size(ndim,1),1);
        tmp = ndim;
        tmp(~ind,:) = tmp(~ind,:) + randn(sum(~ind),2)*0.1;
        s = scatter(tmp(:,1),tmp(:,2));
        lim = [0.5 max(ndim(:))+0.5];
        plot(ax,lim,lim,'k--');
        set(ax,'XLim',lim,'YLim',lim);
        set(s,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',0.3,'SizeData',50); 
        xlabel(ax,'identified dim CSA','Color',Colors(1,:));
        ylabel(ax,{'identified dim','Recurrent CSA'},'Color',Colors(2,:));
        set(gcf,'Position',[0 0 280 240]);
        set(ax,'FontSize',18);
        ef;
    

% fig 2: 1 simulation
    % run simulation to find the good plot
    Magnitudes = [1 3];
    Ranks = [2 2];
    nt = 5; % number of time points in a trial
    nx = 12; ny = 6;
    signse = 1;
    while true
        seed = rng('Shuffle');
        load tmp;rng(seed);
        [X_fit,Y_hist_fit,Y_fit,tr_ind,flts,params] = simu.AR1(struct('rnk',Ranks,'magnitude',Magnitudes,'signse',signse,'ny',ny));
        

        % run CSA
        [B,ops] = regressors.bilinear.train([{X_fit} Y_hist_fit],Y_fit,[2 2]);
        Ypred = regressors.bilinear.predict([{X_fit} Y_hist_fit],B,ops);

        % run RRR
        [B_rrr,ops_rrr] = regressors.bilinear.train([{X_fit} Y_hist_fit],Y_fit,[2 0]);
        Ypred_rrr = regressors.bilinear.predict([{X_fit} Y_hist_fit],B_rrr,ops_rrr);

        % plot predicted by line
            % scale = params.signse / .08;
            % Colors = cbrewer2('Set1',2);
            % ax = np;
            % for ii = 1:size(Y_fit,2)
            %     h = plot(ax,Y_fit(:,ii)+ii*scale,'Color',[0 0 0],'LineWidth',1);
            %     h(2) = plot(ax,Ypred_rrr(:,ii)+ii*scale,'-','Color',Colors(1,:),'LineWidth',.75);
            %     h(3) = plot(ax,Ypred(:,ii)+ii*scale,'-','Color',Colors(2,:),'LineWidth',.75);
            % end
            % set(ax,'XLim',[20 50],'YLim',[0 (size(Y_fit,2)+1)*scale]);
            % l = legend(h,{'true','predicted w/o history','predicted w/ history'},'Location','northoutside','Box','off','FontSize',9);
            % ax.XAxis.Visible = 'off';
            % ax.YAxis.Visible = 'off';
            % set(ax,'FontSize',9);
            % % ax.Position(2) = 0;
            % set(gcf,'Position',[0 0 1.6 1.6]*72);
            % l.Position(1) = 0;
            % l.Position(2) = 1-l.Position(4);
            % ax.Position(2) = 0;
            % ax.Position(4) = l.Position(2) - ax.Position(2);
            % ef;
            % save tmp seed;
        % plot first axis of U / V
            [Us,Ss,Vs] = cellfun(@(x) svd(x,"vector"),[arrayfun(@(s) s.wfilt, flts,'uni',0); B_rrr; B],'uni',0);
            Us = cat(3,Us{:,1});
            Vs = cat(3,Vs{:,1});
            Colors = [0 0 0; cbrewer2('Set1',2)];
            ax = np(2,1);
            % plot true
            % plot(ax(1),flts(1).wx(:,1) / norm(flts(1).wx(:,1)),'.-','LineWidth',1,'Color',[0 0 0]);
            % plot(ax(2),flts(1).wy(1,:) / norm(flts(1).wy(1,:)),'.-','LineWidth',1,'Color',[0 0 0]);
            % plot predicted
            [~,h] = plt_cmp_vectors(squeeze(Us(:,1,:)),[],'lines',ax(1));
            arrayfun(@(ii) set(h(ii),'Marker','.','Color',Colors(ii,:)), 1:numel(h));
            [~,h] = plt_cmp_vectors(squeeze(Vs(:,1,:)),[],'lines',ax(2));
            arrayfun(@(ii) set(h(ii),'Marker','.','Color',Colors(ii,:)), 1:numel(h));

            set(ax,'YLim',[-1 1],'YTick',[-1 0 1],'XTick',[],'FontSize',8);
            % arrayfun(@(ii) set(ax(ii),'XLim',[0 size(Us{1,ii},1)]+0.5), 1:2);
            arrayfun(@(ii) ylabel(ax(ii),'A.U.','FontSize',9), 1:2);
            xlabel(ax(1),'elements in U','FontSize',9);
            xlabel(ax(2),'elements in V','FontSize',9);
            
            set(gcf,'Position',[0 0 1.3 1.6]*72);

            ef;
        pause(5);
    end

    % plot first axis of U / V
    [Us,Ss,Vs] = cellfun(@(x) svd(x,"vector"),[arrayfun(@(s) s.wfilt, flts,'uni',0); B_rrr; B],'uni',0);
    Colors = [0 0 0; cbrewer2('Set1',2)];
    ax = np(2,1);
    % plot true
    % plot(ax(1),flts(1).wx(:,1) / norm(flts(1).wx(:,1)),'.-','LineWidth',1,'Color',[0 0 0]);
    % plot(ax(2),flts(1).wy(1,:) / norm(flts(1).wy(1,:)),'.-','LineWidth',1,'Color',[0 0 0]);
    % plot predicted
    for ii = 1:3
        plot(ax(1),Us{ii,1}(:,1),'.-','LineWidth',1,'Color',Colors(ii,:));
        plot(ax(2),Vs{ii,1}(:,1),'.-','LineWidth',1,'Color',Colors(ii,:));
    end
    set(ax,'YLim',[-1 1],'YTick',[-1 0 1],'XTick',[]);
    arrayfun(@(ii) set(ax(ii),'XLim',[0 size(Us{1,ii},1)]+0.5), 1:2);
    arrayfun(@(ii) ylabel(ax(ii),'A.U.'), 1:2);
    xlabel(ax(1),'elements in U');
    xlabel(ax(2),'elements in V');
    set(gcf,'Position',[0 0 150 200]);

    ef;


% fig 2: simulations, reacreate rank
    % Nsamples = [100 1000 10000 1e5]/4;
    Magnitudes = [1 3];
    Ranks = [2 2];
    nt = 5; % number of time points in a trial
    nx = 12; ny = 10;
    signse = 1;
    Test_ranks = {0:5,[0 5]};
    
    % Nsamples = round((100:100:1000)/nt);
    Nsamples = round(1e3/nt);
    Nsims = 20;
    Errs = nan(Nsims,numel(Test_ranks{1}),numel(Test_ranks{2}),10);
    % load figure_Filters; flts = Filters(14,:);% clear Filters;
    parfor isim = 1:Nsims
        rng('Shuffle');
        flts = gen_flts([nx ny;ny ny],Ranks,Magnitudes);
        % replace wy in history filter; requires having the same rank
        % flts(2).wy(1:Ranks(1),:) = flts(1).wy / sqrt(Magnitudes(1)) * sqrt(Magnitudes(2));
        % flts(2).wfilt = flts(2).wx * flts(2).wy;
        switch 1
        case 1
            [X_fit,Y_hist_fit,Y_fit,tr_ind,~,params] = simu.AR1(struct('rnk',Ranks,...
                                                                    'signse',signse,...
                                                                    'magnitude',Magnitudes,...
                                                                    'ntr',max(Nsamples),'nt',nt,...
                                                                    'nx',size(flts(1).wx,1),'ny',size(flts(1).wy,2)),flts);
        case 2
            % [X_fit,Y_hist_fit,Y_fit,tr_ind,flts,params] = simu.common_input(struct());
        end


        % see the rank of inter-region communication w/ or w/o history
        % run cross validation
        [err,ops_cv,to_save] = myCV.by_trial({X_fit,Y_hist_fit{1}},Y_fit,tr_ind,struct('lambda',0,'rnk',flip(combinations(Test_ranks{2},Test_ranks{1}),2),'n_fold',10));
        Errs(isim,:,:,:) = reshape(err,numel(Test_ranks{1}),numel(Test_ranks{2}),[]);

        % save filters
        Filters(isim,:) = flts;
        fprintf('sim %d done\n',isim);
    end


    rnk = loss.cal_rnk(squeeze(Errs(:,:,1,:)),0:5);
    rnk(end+1,:) = loss.cal_rnk(squeeze(Errs(:,:,2,:)),0:5)
    % err: nsessions / subjects x nranks x nfold
    ax = np;scatter(rnk(1,:)+randn(1,size(rnk,2))*0.1,rnk(2,:)+randn(1,size(rnk,2))*0.1);
    scatter(ax,Ranks(1),Ranks(1),'.');
    equalize_plot(ax);
    xlabel('w/o hist dim');
    ylabel('w/ hist dim');
    ef;
    return

    % plot average error, session by session
    to_plot = find(rnk(1,:)>rnk(2,:))
    n_sessions = size(Errs,1);
    ax = np(n_sessions);
    Colors = cbrewer2('Set1',2);
    % arrayfun(@(ii) plot_multiple_lines(squeeze(1-Errs(ii,:,1,:))',ax(ii),'base_color',Colors(1,:)), 1:n_sessions);
    % arrayfun(@(ii) plot_multiple_lines(squeeze(1-Errs(ii,:,2,:))',ax(ii),'base_color',Colors(2,:)), 1:n_sessions);
    for ii = to_plot
        yyaxis(ax(ii),'left');
        plot_multiple_lines(squeeze(1-Errs(ii,:,1,:))',ax(ii),'base_color',Colors(1,:),'to_plt',{'err'});
        yyaxis(ax(ii),'right');
        plot_multiple_lines(squeeze(1-Errs(ii,:,2,:))',ax(ii),'base_color',Colors(2,:),'to_plt',{'err'});
        ax(ii).YAxis(2).Color = Colors(2,:);
        ax(ii).YAxis(1).Color = Colors(1,:);
    end
    ef;
    

    return
    

    iSession = 15;
    ax = np(2);
    Colors = cbrewer2('Set1',2);
    % yyaxis(ax,'left');
    plot_multiple_lines(squeeze(1-Errs(iSession,:,1,:))',ax(1),'base_color',Colors(1,:),'to_plt',{'err'});
    % yyaxis(ax,'right');
    plot_multiple_lines(squeeze(1-Errs(iSession,:,2,:))',ax(2),'base_color',Colors(2,:),'to_plt',{'err'});
    
    ax.YAxis(2).Color = Colors(2,:);
    ax.YAxis(1).Color = Colors(1,:);
    ef;
    


% fig 2: find one good semedo session
    [f,SessionNames] = load_data.semedo_results('cv_grid_sepstim_semedo');
    Err = cat(1,f(:).Err);

    % average across dataset
    ops = f(1).ops;
    Colors = cbrewer2('Set1',4);


    % session by session
        iSession = 12; % 39 for semedo, 12 for alex data
        ax = np;
        h = plot_multiple_lines(1-squeeze(Err(iSession,2:end,1,:))',ax,'x',ops.cv.Ranks_x(2:end),'base_color',Colors(1,:),'to_plt',{'err'});
        h.err.Marker = 'o';
        h.err.LineWidth = 1.5;
        h = plot_multiple_lines(1-squeeze(Err(iSession,:,end,:))',ax,'x',ops.cv.Ranks_x,'base_color',Colors(2,:),'to_plt',{'err'});
        h.err.Marker = 'o';
        h.err.LineWidth = 1.5;
        % plot_multiple_lines(1-squeeze(Err(iSession,:,end,:))',ax(2),'x',ops.cv.Ranks_x,'base_color',Colors(2,:));
        % title(ax(1),'CSA');
        % title(ax(2),'Recurrent CSA');
        set(ax,'FontSize',9,'XTick',0:5,'XLim',[0 5],'XTickLabelRotation',0);
        xlabel(ax,'input dimensionality','FontSize',10);
        ylabel(ax,'predictability','FontSize',10);
        set(gcf,'Position',[0 0 1.7 2]*72);
        title('example session','FontSize',10);
        ef;



% checked
    % fig 1: schematic simulation
        n_samples = 20;
        rng(8);
        a = normrnd(0,0.5,n_samples,1);
        b = normrnd(0,0.5,n_samples,1);
        c = normrnd(0,0.5,n_samples,1);

        x = 0.5 * a + b;
        y = b + 0.1 * c;
        z = 0.1*x - 0.2 * b + c;
        x = a;
        y = b;
        z = c;
        x = x - mean(x);
        y = y - mean(y);
        z = z - mean(z);
        center = mean([x y z],2);

        % plot
        coeff = [1 1 1 ; -1 0.7 -1;0 3 0];
        sz    = [2 1.3;2 1.5;2 1.5]*180;
        % sz    = [421 258;182 122;182 122];
        % cmap  = {'RdGy','RdYlGn','BrBG'};
        lcolor  = {'#c830cc','#4ea6dc','#4775e7'};
        cmap  = {gray};
        msz = [120 18 18];
        % cmap  = {flip(cmap{1},1), cmap{1}(:,[2 1 3]), cmap{1}(:,[3 2 1])};
        % ltype = {'-','--','-.'};
        ip = 1;
        close all; h = [];


        % comm axis
        % scatter3(x,y,z,[],[0.4 0.4 0.4],'filled'); hold on;
        S = scatter3(x,y,z,msz(ip),[x,y,z]*coeff(ip,:)','filled','MarkerEdgeColor',[0 0 0]); hold on;
        c = colorbar; c.Ticks = []; 
        c.Position = [0.93    0.1100    0.04    0.8150];
        
        % h(end+1) = plt_coeff(coeff(ip,:),center,['k' ltype{ip}],'LineWidth',1.5); % communication axis
        center = mean([x,y,z],1);
        h = plot3([-1 1]*coeff(ip,1)+center(1),[-1 1]*coeff(ip,2)+center(2),[-1 1]*coeff(ip,3)+center(3),'-','Color',lcolor{ip},'LineWidth',2);
        % h(end+1) = plt_coeff(coeff(ip,:),center,'-','Color',lcolor{ip},'LineWidth',2); % communication axis
        ax = gca;
        % set(ax,'XLim',[-1.2686 1.2381],'YLim',[-1.1550 1.2027],'ZLim',[-1.0000 1.0354],...
        % 	'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'View',[-73.4862 31.8924]);

        f = gcf; 
        f.Position([3 4]) = sz(ip,:);
        colormap(cmap{1});

        set(ax,'View',[-72 29],'XTick',-1:1,'YTick',-1:1,'ZTick',-1:1,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
        set(ax,'XLim',[-1.2192 1.8]);
        set(ax,'YLim',[-1 1]*1.3);
        set(ax,'ZLim',[-1.1820  1.1877]);
        set(ax,'XTick',linspace(ax.XLim(1),ax.XLim(2),3),'YTick',linspace(ax.YLim(1),ax.YLim(2),3),'ZTick',linspace(ax.ZLim(1),ax.ZLim(2),3));
        ax.LineWidth = 1.5;
        % set(gca,'XTick',-1:1,'YTick',-1:1,'ZTick',-1:1,'XLim',[-1.2686 1.2381],'YLim',[-1.1550 1.2027],'ZLim',[-1.0000 1.0354],'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'View',[-73.4862 31.8924]);
        % save figure as png
        exportgraphics(gcf,'tmp.png','Resolution',300);

    return
    % fig 2: simulations, better error wise
        % run simulation, saved one simulation results that were used to produce figure in data/cmp_error.mat
            irun = 2; fprintf('%d run\n',irun); rng('Shuffle');
            % Nsamples = [500:100:1000 1e4 1e5]/10;
            % Nsamples = 100*(2.^(0:9)) / 10;
            Nsamples = 100*(2.^(2:9));
            % Nsamples = [100:100:1000 1e4 1e5 1e6]/10;
            Nsims = 10;
            % Errs = nan(Nsims,numel(Nsamples),2,10);
            Errs = nan(Nsims,numel(Nsamples),2);
            err_fun = @(a,b) sum((a-b).^2,'all');
            flts = [];
            for isim = 1:Nsims
                [X_fit,Y_hist_fit,Y_fit,tr_ind,flts,params] = simu.AR1(struct('rnk',[1 5],...
                                                                    'signse',0.1,...
                                                                    'magnitude',[1 5],...
                                                                    'ntr',max(Nsamples),'nt',10,'self_couple',0),flts);
                                                                    % 'ntr',1,'nt',max(Nsamples),'self_couple',0),flts);
                                                                    
                
                for iSampleSize = 1:numel(Nsamples)
                    % original, subsample by trial
                    ind = ismember(tr_ind,1:Nsamples(iSampleSize));
                    % new, subsample by time points
                    % ind = 1:Nsamples(iSampleSize);

                    % estimation error
                    [B,ops_] = bilinear.train({X_fit(ind,:),Y_hist_fit{1}(ind,:)},Y_fit(ind,:),[params.rnk(1) 0]);
                    [B(2,:),ops_] = bilinear.train({X_fit(ind,:),Y_hist_fit{1}(ind,:)},Y_fit(ind,:),params.rnk);
                    Errs(isim,iSampleSize,:) = cellfun(@(x) err_fun(x,flts(1).wfilt), B(:,1));

                    % cv error
                    % [Errs(isim,iSampleSize,:,:),ops,to_save] = myCV.by_trial({X_fit(ind,:),Y_hist_fit{1}(ind,:)},Y_fit(ind,:),tr_ind(ind),struct('rnk',[params.rnk(1) 0; params.rnk]));
                end

                fprintf('sim %d done\n',isim);
            end
            params.seed = rng;
            params.Bmag = sum(flts(1).wfilt.^2,'all');
            save(sprintf('cmp_error_run%d.mat',irun),'Errs','Nsamples','Nsims','params','flts');

            return


        % load and plot
            load('data/cmp_error.mat')
            Colors = [0.8941    0.1020    0.1098; 0.2157    0.4941    0.7216];
            % mean(Errs(:,:,1),1)
            ind = 3:size(Errs,2);
            close all; ax = subplot(1,1,1,'NextPlot','add');
            for jj = 1:2
                h = plot_multiple_lines(Errs(:,ind,jj)/params.Bmag,ax,'x',Nsamples(ind)*params.nt,'base_color',Colors(jj,:),'to_plt',{'err'});
                h.err.LineWidth = 2;
            end
            ax.XAxis.TickLabelRotation = 0;
            xlim([0 max(Nsamples(ind)*params.nt)]);
            xlabel(ax,'# samples');
            ylabel(ax,'Error');
            set(gcf,'Position',[0 0 2.2 1.5]*72);
            set(ax,'FontSize',9);
            

            exportgraphics(gcf,'tmp.png','Resolution',300);

            return


