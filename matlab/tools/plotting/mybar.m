function h = mybar(X,ax)
% bar plot with error bar and individual lines
% X: samples * Nvariables

if nargin < 2
    ax = np;
end

N = size(X,2);

% bar
h.bar    = bar(ax,1:N,mean(X,1),'FaceColor',[0 0 0],'FaceAlpha',0.5,'EdgeColor',[0 0 0]);
% individual sample lines
h.lines  = plot(ax,1:N,X,'-','LineWidth',0.5,'Color',[0 0 0 0.5]);
% errorbar
h.errbar = errorbar(ax,1:N,mean(X,1),std(X,[],1)/sqrt(size(X,1)),'.','Color',[0 0 0],'LineWidth',1);


% figure settings
set(ax,'XTick',1:N);



end