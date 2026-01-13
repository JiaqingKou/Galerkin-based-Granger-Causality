function [X,Y,Mesh2Vec,Vec2Mesh]=CavityGrid(N)
    x=linspace(-1,1,N+1);
    y=x;
    [X,Y]=meshgrid(x,y);
    Mesh2Vec = @(X,Y) [X(:);Y(:)];
    Vec2Mesh = @(xy)  deal(reshape( xy(1:end/2),sqrt(size(xy,1)/2),sqrt(size(xy,1)/2)), ...
                         reshape( xy(end/2+1:end),sqrt(size(xy,1)/2),sqrt(size(xy,1)/2)) );
end