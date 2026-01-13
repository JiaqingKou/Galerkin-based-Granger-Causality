%% Function: visualiztion of the shear flow model (in MATLAB)
%   INPUT
%     Vec9£ºcoefficients of 9 modes in vector, size(1,9)
%
%   OUTPUT
%     plot single 3D fluid snapshot in matlab 
%
function Shearflow3Dflow(Vec9)

    % constants define in modes9.m
    Lx = 4*pi;      % length of the domain in the x direction 
    Lz = 2*pi;      % length of the domain in the z direction
    gridx = 10;     % number of grids in the x direction
    gridy = 20;     % number of grids in the y direction
    gridz = 20;     % number of grids in the z direction
    
    % restruct velocity field base on 9 modes
    [X, Y, Z, u, v, w] = ShearflowVelocity(Vec9, gridx, gridy, gridz);
    % plot velocity field 
    colormap(jet(64));
    % Scale the vectors for better visualization
    scale_factor = 1;
    % Plot the velocity vectors using quiver3
    quiver3(X, Z, Y, u, w, v, scale_factor,'LineWidth', 1, 'color', 'black');
    hold on;
    % Calculate velocity magnitude
    velocity_magnitude = sqrt(u.^2 + v.^2 + w.^2);
    % Create a 3D scatter plot with color-coded points based on velocity magnitude
    scatter3(X(:), Z(:), Y(:), 10, velocity_magnitude(:), 'filled');
    % Plot the velocity magnitude on the surfaces of the cube
    surface_color = squeeze(velocity_magnitude(:, 1, :));
    h = surf(squeeze(X(:, end, :)), squeeze(Z(:, end, :)), squeeze(Y(:, end, :)), surface_color);
    set(h, 'EdgeColor', 'interp', 'FaceColor', 'interp');
    surface_color = squeeze(velocity_magnitude(1, :, :));
    h = surf(squeeze(X(1, :, :)), squeeze(Z(1, :, :)), squeeze(Y(1, :, :)), surface_color);
    set(h, 'EdgeColor', 'interp', 'FaceColor', 'interp');
    surface_color = squeeze(velocity_magnitude(:, :, 1));
    h = surf(squeeze(X(:, :, end)), squeeze(Z(:, :, end)), squeeze(Y(:, :, end)), surface_color);
    set(h, 'EdgeColor', 'interp', 'FaceColor', 'interp');
    % Create edges of the cube
    p = [0 -1 0; 0 -1 Lz; 0 1 Lz; 0 1 0; 0 -1 0;...
         Lx -1 0; Lx -1 Lz; Lx 1 Lz; Lx 1 0; Lx -1 0;...
         Lx 1 0; 0 1 0; 0 1 Lz; Lx 1 Lz; Lx -1 Lz; 0 -1 Lz];
    x = p(:,1)';
    y = p(:,2)';
    z = p(:,3)';
    for i=1:15
        line(x(i:i+1),z(i:i+1),y(i:i+1),'Color','black','LineWidth',3);
    end
    line([0 Lx],[Lz/2 Lz/2],[0 0],'Color','black','LineWidth',2,'LineStyle','--');
    line([Lx/2 Lx/2],[Lz/2 Lz/2],[-1 1],'Color','black','LineWidth',2,'LineStyle','--');
    line([Lx/2 Lx/2],[0 Lz],[0 0],'Color','black','LineWidth',2,'LineStyle','--');

    hold off;
    xlabel('X-axis', 'FontSize',12, 'FontAngle','italic');
    ylabel('Z-axis', 'FontSize',12, 'FontAngle','italic');
    zlabel('Y-axis', 'FontSize',12, 'FontAngle','italic');
    title('3D Fluid Contour', 'FontSize',14, 'FontWeight','bold');
    % Adjust the size of the fig window
    width = 1000;
    height = 600;
    set(gcf, 'Position', [100, 0, width, height]);
    zlim([-2,2]);

    % Add colorbar
    colorbar;
    ylim([0,Lz+1])
    
    % plot 9-mode factors
    labels = {'Mode1','Mode2','Mode3';'Mode4','Mode5','Mode6';'Mode7','Mode8','Mode9'};
    for i = 1:3
        for j = 1:3
            data = reshape(Vec9,[3,3])';
            if i==1 && j==1 && data(1,1)>0.3
                text('String', num2str(data(i,j),'%.2f'), 'Position', [-0.2+0.075*i 1.0-0.075*j], ...
                     'Units','normalized', 'HorizontalAlignment','left', 'FontSize',15, 'color','red',...
                     'VerticalAlignment','bottom', 'FontAngle','italic', 'FontWeight','bold');
                text('String', labels{i,j}, 'Position', [-0.2+0.075*i 1.04-0.075*j], ...
                     'Units','normalized', 'HorizontalAlignment','left', 'FontSize',10, ...
                     'VerticalAlignment','bottom',  'FontWeight','bold', 'color','red');
                continue
            end
            text('String', num2str(data(i,j),'%.2f'), 'Position', [-0.2+0.075*i 1.0-0.075*j], ...
                 'Units','normalized', 'HorizontalAlignment','left', 'FontSize',15, ...
                 'VerticalAlignment','bottom', 'FontAngle','italic', 'FontWeight','bold');
            text('String', labels{i,j}, 'Position', [-0.2+0.075*i 1.04-0.075*j], ...
                 'Units','normalized', 'HorizontalAlignment','left', 'FontSize',10, ...
                 'VerticalAlignment','bottom',  'FontWeight','bold');
        end
    end
    annotation('textbox', [0.02, 0.885, 0.1, 0.1], 'String', 'Mode Coefficient.:', ...
               'EdgeColor', 'none', 'FontSize', 16, 'FontWeight', 'bold', ...
               'color','red', 'FontAngle','italic');

end