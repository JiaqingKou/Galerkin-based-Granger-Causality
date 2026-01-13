%%  Causality-guided flow control to suppress extreme events
%   DISCRIPTION:
%     upon detection of an extreme event (a_1 > a_th), 
%     a control body force is applied to the modulator (mode 3) to facilitate energy transfer
%     from mode 1 to mode 2, thereby suppressing the extreme event.
%
%   REFERENCE:
%     see §5.4 of [Wang, Y., Kou, J., Noack, B.R. et al. (2026). Spectral analysis of mixing in 2D high-Reynolds ﬂows]
%     and [Moehlis, J., Faisst, H., & Eckhardt, B. (2004). A low-dimensional model for turbulent shear flows]
%
clear all
addpath('./Code_PlotResult');

% model parameters 
Re = 600;       % Reynolds number
Lx = 4*pi;      % length of the domain in the 𝑥 direction 
Lz = 2*pi;      % length of the domain in the 𝑧 direction
Alpha = 2*pi/Lx;
Beta = pi/2;
Gamma = 2*pi/Lz;
K_ag = sqrt(Alpha^2 + Gamma^2);
K_bg = sqrt(Beta^2 + Gamma^2);
K_abg = sqrt(Alpha^2 + Beta^2 + Gamma^2);
Time = 5000;    % total time of the control
      
% initial condition
IC = [0; 0.03; 0.1; 0.05; 0.1; 0; 0; 0.04; 0];

%% Causality-guide control on suppression of extreme events
% control parameters 
K = 0.01;       % control amplitude
a_th = 0.4;     % control objective value for mode 1 / activation value for control
dtc = 0.1;      % time step of flow control    

% perform flow control
disp('simulating flow control...')
A = zeros(Time/dtc+1,9);    % time series of mode coefficients 
A(1,:) = IC;
pairmode = 2;               % modulated pair mode with mode 1
fc = 0;                     % applied control body force
for i = 1:Time/dtc
    % the control force fc is injected into the dynamic equation of a(3)
    fun = @(t,a)  [Beta^2/Re*(1-a(1)) - sqrt(1.5)*Beta*Gamma/K_abg*a(6)*a(8) + sqrt(1.5)*Beta*Gamma/K_bg*a(2)*a(3);....          
                  -(4/3*Beta^2+Gamma^2)/Re*a(2) + 5/3*sqrt(2/3)*Gamma^2/K_ag*a(4)*a(6) - Gamma^2/sqrt(6)/K_ag*a(5)*a(7) - Alpha*Beta*Gamma/sqrt(6)/K_ag/K_abg*a(5)*a(8) - sqrt(1.5)*Beta*Gamma/K_bg*(a(1)*a(3)+a(3)*a(9));...
                  fc-(Beta^2+Gamma^2)/Re*a(3) + 2/sqrt(6)*Alpha*Beta*Gamma/K_ag/K_bg*(a(4)*a(7)+a(5)*a(6)) + (Beta^2*(3*Alpha^2+Gamma^2)-3*Gamma^2*(Alpha^2+Gamma^2))/sqrt(6)/K_ag/K_bg/K_abg*a(4)*a(8);...
                  -(3*Alpha^2+4*Beta^2)/3/Re*a(4) - Alpha/sqrt(6)*a(1)*a(5) - 10/3/sqrt(6)*Alpha^2/K_ag*a(2)*a(6) - sqrt(1.5)*Alpha*Beta*Gamma/K_ag/K_bg*a(3)*a(7) - sqrt(1.5)*Alpha^2*Beta^2/K_ag/K_bg/K_abg*a(3)*a(8) - Alpha/sqrt(6)*a(5)*a(9);...
                  -(Alpha^2+Beta^2)/Re*a(5) + Alpha/sqrt(6)*(a(1)*a(4)+a(4)*a(9)) + Alpha^2/sqrt(6)/K_ag*a(2)*a(7) - Alpha*Beta*Gamma/sqrt(6)/K_ag/K_abg*a(2)*a(8) + 2/sqrt(6)*Alpha*Beta*Gamma/K_ag/K_bg*a(3)*a(6);...
                  -(3*Alpha^2+4*Beta^2+3*Gamma^2)/3/Re*a(6) + Alpha/sqrt(6)*a(1)*a(7) + sqrt(1.5)*Beta*Gamma/K_abg*a(1)*a(8) + 10/3/sqrt(6)*(Alpha^2-Gamma^2)/K_ag*a(2)*a(4) - 2*sqrt(2/3)*Alpha*Beta*Gamma/K_ag/K_bg*a(3)*a(5) + Alpha/sqrt(6)*a(7)*a(9) + sqrt(1.5)*Beta*Gamma/K_abg*a(8)*a(9);...
                  -(Alpha^2+Beta^2+Gamma^2)/Re*a(7) - Alpha/sqrt(6)*(a(1)*a(6)+a(6)*a(9)) + 1/sqrt(6)*(Gamma^2-Alpha^2)/K_ag*a(2)*a(5) + 1/sqrt(6)*Alpha*Beta*Gamma/K_ag/K_bg*a(3)*a(4);...
                  -(Alpha^2+Beta^2+Gamma^2)/Re*a(8) + 2/sqrt(6)*Alpha*Beta*Gamma/K_ag/K_abg*a(2)*a(5) + Gamma^2*(3*Alpha^2-Beta^2+3*Gamma^2)/sqrt(6)/K_ag/K_bg/K_abg*a(3)*a(4);...
                  -9*Beta^2/Re*a(9) + sqrt(1.5)*Beta*Gamma/K_bg*a(2)*a(3) - sqrt(1.5)*Beta*Gamma/K_abg*a(6)*a(8)];
    % solve equations for one time step
    [~,temp] = ode45(fun, [0 dtc], A(i,:)');
    A(i+1,:) = temp(end,:);
    a_1 = A(i+1,1);
    % apply control when extreme event occurs
    if a_1 > a_th
        fc = -sign(A(i+1,pairmode)) * K;    % bang-bang control
    else
        fc = 0;                             % no control                       
    end
    % progress visualization
    if mod(i*dtc,1000) == 0
       disp(['step ',num2str(i),'/',num2str(Time/dtc),' simulated...']) 
    end
end

%% Results Visualization 
% % plot mode coefficients (in Matlab) 
figure('Position', [100, 100, 1200, 400]); 
plot([0:dtc:Time],A);
line([0 Time], [a_th a_th], 'Color','red','LineStyle','--')
xlabel('time', 'FontSize',16,'FontWeight','bold');
ylabel('mode coefficient', 'FontSize',16,'FontWeight','bold');
title('Evolution of 9-mode shear flow', 'FontSize',20,'FontWeight','bold')
grid on;
legend('a1','a2','a3','a4','a5','a6','a7','a8','a9');

% % plot flow contours (in Matlab) 
figure; Shearflow3Dflow(A(1,:))

% % plot flow contours (in Tecplot) 
% ShearflowTecplot(A(1,:),dtc,'tecplot')



