function X = vec2mats(v,s1,s2)
% transforming one vector into multiple matrices
% Input: v: vector
%        s1: size of the first dimension of the matrices
%        s2: size of the second dimension of the matrices

nmat = numel(s1);
X = cell(1,nmat);
ind_start = 1;
for ii = 1:nmat
    n_elements = s1(ii)*s2(ii);
    X{ii} = reshape(v(ind_start:(ind_start+n_elements-1)),s1(ii),s2(ii));
    ind_start = ind_start + n_elements;
end


end