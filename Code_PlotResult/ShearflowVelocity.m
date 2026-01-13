%% Function: restruct velocity based on Galerkin modes
%   INPUT
%     vec9£ºrow array of 9-mode coefficients, size(1,9)
%     gridx/y/z: number of grids in x,y,z direction, respectively
%
%   OUTPUT
%     x/y/z_mat: x,y,z coordinates on all grid points, size(nx,ny,nz)
%     u/v/w_mat: velocity components on all grid points, size(nx,ny,nz)
%
%   REFERENCE:
%     see [Moehlis, J., Faisst, H., & Eckhardt, B. (2004). A low-dimensional model for turbulent shear flows]
%     and [Wang, Y., Kou, J., Noack, B.R. et al. (2026). Causal analysis of a turbulent shear flow model]
%
%% function
function [x_mat,y_mat,z_mat,u_mat,v_mat,w_mat] = ShearflowVelocity(vec9, gridx, gridy, gridz)
    % optional inputs
    if isempty(gridx) gridx=40; end
    if isempty(gridy) gridy=20; end
    if isempty(gridz) gridz=10; end
    % constants of the shear flow model
    Lx = 4*pi;      % length of the domain in the x direction 
    Lz = 2*pi;      % length of the domain in the z direction
    Alpha = 2*pi/Lx;
    Beta = pi/2;
    Gamma = 2*pi/Lz;
    x = linspace(0, Lx, gridy+1);   % x-coordinates
    y = linspace(-1, 1, gridx+1);   % y-coordinates
    z = linspace(0, Lz, gridz+1);   % z-coordinates
    [x_mat, y_mat, z_mat] = meshgrid(x, y, z);  % Create a denser 3D grid
    u_mat = zeros(size(x_mat));
    v_mat = zeros(size(x_mat));
    w_mat = zeros(size(x_mat));
    if size(vec9) ~= [1,9]
        error('function velocity(modes): input must be an array of size(1,9)');
    end
    % calculation of Fourier modes (see appendix B of the reference)
    if vec9(1)
        u_mat = u_mat + vec9(1)*sqrt(2).*sin(pi*y_mat/2);
    end
    if vec9(2)
        u_mat = u_mat + vec9(2)*4/sqrt(3).*cos(pi*y_mat/2).*cos(pi*y_mat/2).*cos(Gamma*z_mat);
    end
    if vec9(3)
        temp = 2/sqrt(4*Gamma^2+pi^2);
        v_mat = v_mat + vec9(3)*temp*2*Gamma.*cos(pi*y_mat/2).*cos(Gamma*z_mat);
        w_mat = w_mat + vec9(3)*temp*pi.*sin(pi*y_mat/2).*sin(Gamma*z_mat);
    end
    if vec9(4)
        w_mat = w_mat + vec9(4)*4/sqrt(3).*cos(Alpha*x_mat).*cos(pi*y_mat/2).*cos(pi*y_mat/2);
    end
    if vec9(5)
        w_mat = w_mat + vec9(5)*2.*sin(Alpha*x_mat).*sin(pi*y_mat/2);
    end
    if vec9(6)
        temp = 4*sqrt(2)/sqrt(3*(Alpha^2+Gamma^2));
        u_mat = u_mat + vec9(6)*temp*(-Gamma).*cos(Alpha*x_mat).*cos(pi*y_mat/2).*cos(pi*y_mat/2).*sin(Gamma*z_mat);
        w_mat = w_mat + vec9(6)*temp*Alpha.*sin(Alpha*x_mat).*cos(pi*y_mat/2).*cos(pi*y_mat/2).*cos(Gamma*z_mat);
    end
    if vec9(7)
        temp = 2*sqrt(2)/sqrt(Alpha^2+Gamma^2);
        u_mat = u_mat + vec9(7)*temp*Gamma.*sin(Alpha*x_mat).*sin(pi*y_mat/2).*sin(Gamma*z_mat);
        w_mat = w_mat + vec9(7)*temp*Alpha.*cos(Alpha*x_mat).*sin(pi*y_mat/2).*cos(Gamma*z_mat);
    end
    if vec9(8)
        N8 = 2*sqrt(2)/sqrt((Alpha^2+Gamma^2)*(4*Alpha^2+4*Gamma^2+pi^2));
        u_mat = u_mat + vec9(8)*N8*pi*Alpha.*sin(Alpha*x_mat).*sin(pi*y_mat/2).*sin(Gamma*z_mat);
        v_mat = v_mat + vec9(8)*N8*2*(Alpha^2+Gamma^2).*cos(Alpha*x_mat).*cos(pi*y_mat/2).*sin(Gamma*z_mat);
        w_mat = w_mat + vec9(8)*N8*(-pi)*Gamma.*cos(Alpha*x_mat).*sin(pi*y_mat/2).*cos(Gamma*z_mat);
    end
    if vec9(9)
        u_mat = u_mat + vec9(9)*sqrt(2).*sin(3*pi*y_mat/2);
    end
end