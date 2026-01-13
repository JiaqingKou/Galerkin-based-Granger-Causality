function [ h,mycolormap ] = PlotVorticity( STATE,labels,varargin )
    mycolormap = VorticityColormap( 64 );
            h=0;
    if nargout>1
        return;
    end
    % cavity grid and PDE operators
    N = sqrt(size(STATE,1))-1;

    [Grid]=CollocationGrid_q(N);    % the computational grid
    [~,Grid.wc] = clencurt(N);Grid.W = kron(Grid.wc,Grid.wc);   % integration weights
     disp('grid created')

    [ Operators ] = CreateOperators_psi( Grid.D,eye(N+1));
    disp('operators for psi created')

    % grid for the plot
    [xx2,yy2] = meshgrid(Grid.x,Grid.x);
    [xxx,yyy] = meshgrid(-1:.005:1,-1:.005:1); 
    % vorticity 
    [~,frame] = size(STATE);
    w = - Operators.del2*STATE;
    ww = reshape(real(w),N+1,N+1,frame);
    for i=1:frame
        www(:,:,i) = interp2(xx2,yy2,ww(:,:,i),xxx,yyy,'spline');
    end

    a = 10;
    if ~isempty(varargin)
        a = varargin{1};
    end

    www(www>a)=a;www(www<-a)=-a;
    if frame > 3
        row = 2;
    else
        row = 1;
    end
    figure;
    for i=1:frame
        subplot(row,round(frame/row),i)
        contourf(xxx,yyy,www(:,:,i),100,'LineStyle','None') ;    % vorticty
        axis square; 
        box on;
        colormap(mycolormap);
        caxis([-a,a]);
        ax = gca;
        ax.XTick = [-1  1];
        ax.YTick = [-1  1];
        title(['POD mode NO.',num2str(i)], 'FontWeight','bold','FontSize',24);
    end
end

function [ cmap ] = VorticityColormap( ncol )
    %VORTICITYCOLORMAP Summary of this function goes here
    %   Detailed explanation goes here
    cmap=jet(ncol);
    for i=ncol/8+1:7*ncol/16
        cmap(i,:)=[(i-ncol/8)/(5*ncol/16),(i-ncol/8)/(5*ncol/16),1];    % blue fading to white
    end
        cmap(7*ncol/16+1:9*ncol/16,:)=1;  % white
    for i=9*ncol/16+1:7*ncol/8
        cmap(i,:)=[1,1-(i-9*ncol/16)/(5*ncol/16),1-(i-9*ncol/16)/(5*ncol/16)];    % white turning red
    end
end