classdef bilinear < handle
% bilinear regressor



	% train and predict model
	methods (Static)

		function B_full = refill_B(B_part,ops)
			% fill back zeros for exlcuded predictors in B
			% usually excluded by setting rnk to 0
			for ii = 1:numel(B_part)
				ind = find(ops.exclude_id{ii});
				B_full{ii} = insertrows(B_part{ii},zeros(1,size(B_part{ii},2)),ind-1);
			end
		end

		function [Y,Y_pred_parts] = predict(X,B,ops)
			% given input matrices X and coefficients B, predict Y
			% Input: X is a cell array of input matrices
			%        B is a cell array of coefficient matrices

			% need to exclude cells
			if nargin > 2 && isfield(ops,'exclude_id')
				X = arrayfun(@(ii) X{ii}(:,~ops.exclude_id{ii}),1:numel(X),'uni',0);
			end

			% center X based on training data
			if nargin > 2 && isfield(ops,'Xmean')
				X = arrayfun(@(ii) X{ii} - ops.Xmean{ii}(:,~ops.exclude_id{ii}),1:numel(X),'uni',0);
			end

			% calculate predicted Y
			Y = zeros(size(X{1},1),size(B{1},2));
			Y_pred_parts = nan(length(X),size(Y,1),size(Y,2));
			for ii = 1:length(X)
				Y_pred_parts(ii,:,:) = X{ii}*B{ii};
				Y = Y + squeeze(Y_pred_parts(ii,:,:));
			end
			Y = Y + ops.Ymean;

		end

		function [B,ops] = train(X,Y,rnk,ops)
			% train bilinear regressor
			% Input: X is a cell array of input matrices
			%        Y is a matrix of output
			%        rnk is the rank of the coefficient matrices
			% 		 ops is a struct of additional parameters

			if nargin < 4
				ops = struct();
			end

			% demean X
			if ~isfield(ops,'Xmean')
				for ii = 1:numel(X)
					ops.Xmean{ii} = mean(X{ii},1);
				end
			end

			% for each regressor, remove X with zero variance
			for ii = 1:numel(X)
				tmp = var(X{ii},[],1);
				ops.exclude_id{ii} = (tmp==0);

				% remove X with zero variance
				X{ii} = X{ii}(:,~ops.exclude_id{ii}); 
				X{ii} = X{ii} - ops.Xmean{ii}(:,~ops.exclude_id{ii});
			end
			
			% demean y
			ops.Ymean = getOr(ops,'Ymean',mean(Y,1));
			Y = Y - ops.Ymean;

			% check rank is lower than n_cells after removing 0 variance
			for ii = 1:numel(X)
				if size(X{ii},2) < rnk(ii)
					rnk(ii) = size(X{ii},2);
				end
			end



			% use svd method if only one input
			ops.lambda = getOr(ops,'lambda',0);
			if numel(X) == 1 %|| 
				[B{1},wtfit{1},wxfit{1}] = svd_RRR(X{1}, Y, rnk, ops.lambda);
				ops.niter = 1;

			% multiple inputs, but only one has non-zero rnk
			elseif sum(rnk>0)<=1
				
				% find the input with non-zero rnk
				ii = find(rnk>0);
				% run RRR
				if ~isempty(ii)
					[B{ii},wtfit{ii},wxfit{ii}] = svd_RRR(X{ii}, Y, rnk(ii), ops.lambda);
				end

				% set zero rank B as 0
				for ii = find(rnk==0)
					B{ii} = zeros(size(X{ii},2),size(Y,2));
					wtfit{ii} = zeros(size(X{ii},2),1);
					wxfit{ii} = zeros(1,size(Y,2));
				end

				% for extra saving
				ops.niter = 1;


			% multiple inputs
			else
				% easier to just remove rank 0 input here
				% set zero rank B as 0
				for ii = find(rnk==0)
					B{ii} = zeros(size(X{ii},2),size(Y,2));
					wtfit{ii} = zeros(size(X{ii},2),1);
					wxfit{ii} = zeros(1,size(Y,2));
				end

				tic;
				[wtfit(rnk>0),wxfit(rnk>0),wwfilts(rnk>0),to_save] = bilinearMultifiltRRR_coordAscent(X(rnk>0),Y,rnk(rnk>0),ops.lambda,struct('Display','none','MaxIter',1e3));
				ops.niter = to_save.niter;
				tmp = toc;
				% fprintf('Elapsed time: %.2f s after %d iteration\n',tmp, to_save.niter);

				% calculate B
				for ii = find(rnk>0)
					B{ii} = wtfit{ii}*wxfit{ii};
				end
			end

			% refill B?
			if isfield(ops,'refill') && ops.refill
				B = regressors.bilinear.refill_B(B,ops);
			end

			% save other outputs
			ops.wtfit = wtfit;
			ops.wxfit = wxfit;
			ops.rnk   = rnk;

		end

	
	end


end