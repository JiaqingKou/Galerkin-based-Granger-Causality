%% --- function: restructing velocity base on 9 modes ---
% INPUT
%   modes£ºrow array of 9-mode coefficients, size(1,9)
%   x/y/z_mat: x/y/z position of all data points, size(nx,ny,nz)
% OUTPUT
%   u/v/w_mat: u/v/w velocity component of all data points, size(nx,ny,nz)
%
function [u_mat,v_mat,w_mat] = velocity(mode,x_mat,y_mat,z_mat)

    global Lx
    global Lz
    global Alpha
    global Beta
    global Gamma 
    u_mat = zeros(size(x_mat));
    v_mat = zeros(size(x_mat));
    w_mat = zeros(size(x_mat));
    
    if size(mode) ~= [1,9]
        error('function velocity(modes): input must be an array of size(1,9)');
    end
    
    if mode(1)
        u_mat = u_mat + mode(1)*sqrt(2).*sin(pi*y_mat/2);
    end
    if mode(2)
        u_mat = u_mat + mode(2)*4/sqrt(3).*cos(pi*y_mat/2).*cos(pi*y_mat/2).*cos(Gamma*z_mat);
    end
    if mode(3)
        temp = 2/sqrt(4*Gamma^2+pi^2);
        v_mat = v_mat + mode(3)*temp*2*Gamma.*cos(pi*y_mat/2).*cos(Gamma*z_mat);
        w_mat = w_mat + mode(3)*temp*pi.*sin(pi*y_mat/2).*sin(Gamma*z_mat);
    end
    if mode(4)
        w_mat = w_mat + mode(4)*4/sqrt(3).*cos(Alpha*x_mat).*cos(pi*y_mat/2).*cos(pi*y_mat/2);
    end
    if mode(5)
        w_mat = w_mat + mode(5)*2.*sin(Alpha*x_mat).*sin(pi*y_mat/2);
    end
    if mode(6)
        temp = 4*sqrt(2)/sqrt(3*(Alpha^2+Gamma^2));
        u_mat = u_mat + mode(6)*temp*(-Gamma).*cos(Alpha*x_mat).*cos(pi*y_mat/2).*cos(pi*y_mat/2).*sin(Gamma*z_mat);
        w_mat = w_mat + mode(6)*temp*Alpha.*sin(Alpha*x_mat).*cos(pi*y_mat/2).*cos(pi*y_mat/2).*cos(Gamma*z_mat);
    end
    if mode(7)
        temp = 2*sqrt(2)/sqrt(Alpha^2+Gamma^2);
        u_mat = u_mat + mode(7)*temp*Gamma.*sin(Alpha*x_mat).*sin(pi*y_mat/2).*sin(Gamma*z_mat);
        w_mat = w_mat + mode(7)*temp*Alpha.*cos(Alpha*x_mat).*sin(pi*y_mat/2).*cos(Gamma*z_mat);
    end
    if mode(8)
        N8 = 2*sqrt(2)/sqrt((Alpha^2+Gamma^2)*(4*Alpha^2+4*Gamma^2+pi^2));
        u_mat = u_mat + mode(8)*N8*pi*Alpha.*sin(Alpha*x_mat).*sin(pi*y_mat/2).*sin(Gamma*z_mat);
        v_mat = v_mat + mode(8)*N8*2*(Alpha^2+Gamma^2).*cos(Alpha*x_mat).*cos(pi*y_mat/2).*sin(Gamma*z_mat);
        w_mat = w_mat + mode(8)*N8*(-pi)*Gamma.*cos(Alpha*x_mat).*sin(pi*y_mat/2).*cos(Gamma*z_mat);
    end
    if mode(9)
        u_mat = u_mat + mode(9)*sqrt(2).*sin(3*pi*y_mat/2);
    end
end