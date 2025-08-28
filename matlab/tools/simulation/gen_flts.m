function flts = gen_flts(sz,rnk,magnitude)
% generate filters based on size and rank 
% Input: sz -  n filter * 2, size of the filter
%        rnk - n filter * 1, rank of the filter

% initiation
if nargin < 3
    magnitude = ones(size(rnk));
end

% create filter
for ii = 1:numel(rnk)
    % sanity check, rnk should be smaller than the size
    if rnk(ii) > min(sz(ii,:))
        error('rank should be smaller than the size');
    end
    w = randn(sz(ii,1),rnk(ii)) * randn(rnk(ii),sz(ii,2));
    [u,~,v] = svd(w);
    flts(ii).wx = u(:,1:rnk(ii)) * sqrt(magnitude(ii)); % input 
    flts(ii).wy = v(:,1:rnk(ii))' * sqrt(magnitude(ii)); % output
    flts(ii).wfilt = flts(ii).wx * flts(ii).wy; % filter
end

end

% test
% flts = gen_flts([3 5; 5 5],[1 2]);
% ax = np; imagesc(flts.wfilt);colorbar;ef;
