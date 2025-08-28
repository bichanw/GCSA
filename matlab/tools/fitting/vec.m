function v = vec(x,x2)
% VEC - vectorizes the input
% v = vec(x)
%

if nargin == 1
    v = x(:);
elseif nargin == 2
    % wt (x), wx (x2) concatenation
    % both would be cell arrays
    v = catcell(arrayfun(@(ii) reshape(x{ii}*x2{ii},[],1),1:numel(x),'UniformOutput',false),1);
end