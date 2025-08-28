function ard_params = params2ard_true(flts,rnks)
% convert parmas strcut to ard_params struct\
% this probably needs fixing as well

% perform svd on true filters
for jj = 1:numel(flts)
    [u0,s{jj},v0] = svd(flts(jj).wfilt,'econ','vector'); % SVD
    ii = 1:rnks(jj); % indices of singular vectors to keep
    nx_sqrt = sqrt(size(u0,1));
    wU{jj} = u0(:,ii) / nx_sqrt;  % left singular vectors
    wVt{jj} = nx_sqrt * diag(s{jj}(ii))*v0(:,ii)'; % right singular vectors (weighted by singular values)
    s{jj} = s{jj}(ii) * nx_sqrt;
end

ard_params.U_vi_cell = wU';
ard_params.V_vi_cell = wVt';
% ard_params.s = arrayfun(@(i_flt) s{i_flt}(1:rnks(i_flt)) ,1:numel(s),'uni',0);
ard_params.s = s;

% % see what the scales are
% T = 1e3;
% nx = 25;
% U = randn(T,nx) / sqrt(nx);
% ax = np; histogram(vecnorm(U'));ef;

end