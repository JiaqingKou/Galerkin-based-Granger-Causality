%% Function: plot causal map
function causalmap(P,cmin,cmax,cm)

if nargin < 4 || isempty(cm), cm =flipud(bone); end;

n = size(P,1);
assert(ismatrix(P) && size(P,2) == n,'input must be a square 2D matrix');

colormap(cm);
if nargin < 3 || isempty(cmin) || isempty(cmax)
    cmax = max(P(:));
    cmin = min(P(:));
end;
imagesc(P,[cmin cmax]);
axis('square');
xlabel('cause','Fontsize',24,'Fontweight','bold','FontName','Arial');
ylabel('effect','Fontsize',24,'Fontweight','bold','FontName','Arial');
set(gca,'XTick',1:n,'Fontsize',15,'FontName','Arial');
set(gca,'XTickLabel',1:n,'Fontsize',15,'FontName','Arial');
set(gca,'YTick',1:n,'Fontsize',15,'FontName','Arial');
set(gca,'YTickLabel',1:n,'Fontsize',15,'FontName','Arial');
% colorbar

