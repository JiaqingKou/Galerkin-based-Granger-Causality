%%  Galerkin-based Granger Causality (GGC) analysis on turbulent shear flow model
%   DISCRIPTION:
%     Perform causal analysis on the empirical Galerkin model for shear flows 
%     developed by Moehlis et al. (2004), utilizing proposed GGC method.
%     See [Wang, Y., Kou, J., Noack, B.R. et al. (2026). Causal analysis of a turbulent shear flow model]
%
%   DATA SOURCE:
%     [Moehlis, J., Faisst, H., & Eckhardt, B. (2004). A low-dimensional model for turbulent shear flows]
%
clear all
addpath('./Code_PlotResult');

%% --- Turbulent shear flow model ---
% model parameters
Re = 600;       % Reynolds number
Lx = 4*pi;      % length of the domain in x direction 
Lz = 2*pi;      % length of the domain in z direction
Alpha = 2*pi/Lx;
Beta = pi/2;
Gamma = 2*pi/Lz;
K_ag = sqrt(Alpha^2 + Gamma^2);
K_bg = sqrt(Beta^2 + Gamma^2);
K_abg = sqrt(Alpha^2 + Beta^2 + Gamma^2);
Time = 12000;
Dt = 0.1;       % time step of numerical integration
      
% Galerkin-projected 9-mode shear flow ODEs (Eqs. 5.1–5.9) 
f=@(t,a) [Beta^2/Re*(1-a(1)) - sqrt(1.5)*Beta*Gamma/K_abg*a(6)*a(8) + sqrt(1.5)*Beta*Gamma/K_bg*a(2)*a(3);....          
          -(4/3*Beta^2+Gamma^2)/Re*a(2) + 5/3*sqrt(2/3)*Gamma^2/K_ag*a(4)*a(6) - Gamma^2/sqrt(6)/K_ag*a(5)*a(7) - Alpha*Beta*Gamma/sqrt(6)/K_ag/K_abg*a(5)*a(8) - sqrt(1.5)*Beta*Gamma/K_bg*(a(1)*a(3)+a(3)*a(9));...
          -(Beta^2+Gamma^2)/Re*a(3) + 2/sqrt(6)*Alpha*Beta*Gamma/K_ag/K_bg*(a(4)*a(7)+a(5)*a(6)) + (Beta^2*(3*Alpha^2+Gamma^2)-3*Gamma^2*(Alpha^2+Gamma^2))/sqrt(6)/K_ag/K_bg/K_abg*a(4)*a(8);...
          -(3*Alpha^2+4*Beta^2)/3/Re*a(4) - Alpha/sqrt(6)*a(1)*a(5) - 10/3/sqrt(6)*Alpha^2/K_ag*a(2)*a(6) - sqrt(1.5)*Alpha*Beta*Gamma/K_ag/K_bg*a(3)*a(7) - sqrt(1.5)*Alpha^2*Beta^2/K_ag/K_bg/K_abg*a(3)*a(8) - Alpha/sqrt(6)*a(5)*a(9);...
          -(Alpha^2+Beta^2)/Re*a(5) + Alpha/sqrt(6)*(a(1)*a(4)+a(4)*a(9)) + Alpha^2/sqrt(6)/K_ag*a(2)*a(7) - Alpha*Beta*Gamma/sqrt(6)/K_ag/K_abg*a(2)*a(8) + 2/sqrt(6)*Alpha*Beta*Gamma/K_ag/K_bg*a(3)*a(6);...
          -(3*Alpha^2+4*Beta^2+3*Gamma^2)/3/Re*a(6) + Alpha/sqrt(6)*a(1)*a(7) + sqrt(1.5)*Beta*Gamma/K_abg*a(1)*a(8) + 10/3/sqrt(6)*(Alpha^2-Gamma^2)/K_ag*a(2)*a(4) - 2*sqrt(2/3)*Alpha*Beta*Gamma/K_ag/K_bg*a(3)*a(5) + Alpha/sqrt(6)*a(7)*a(9) + sqrt(1.5)*Beta*Gamma/K_abg*a(8)*a(9);...
          -(Alpha^2+Beta^2+Gamma^2)/Re*a(7) - Alpha/sqrt(6)*(a(1)*a(6)+a(6)*a(9)) + 1/sqrt(6)*(Gamma^2-Alpha^2)/K_ag*a(2)*a(5) + 1/sqrt(6)*Alpha*Beta*Gamma/K_ag/K_bg*a(3)*a(4);...
          -(Alpha^2+Beta^2+Gamma^2)/Re*a(8) + 2/sqrt(6)*Alpha*Beta*Gamma/K_ag/K_abg*a(2)*a(5) + Gamma^2*(3*Alpha^2-Beta^2+3*Gamma^2)/sqrt(6)/K_ag/K_bg/K_abg*a(3)*a(4);...
          -9*Beta^2/Re*a(9) + sqrt(1.5)*Beta*Gamma/K_bg*a(2)*a(3) - sqrt(1.5)*Beta*Gamma/K_abg*a(6)*a(8)];
      
% initial condition 
IC = [0; 0.1; 0.1; 0.1; 0.1; 0; 0; 0; 0];

% solve equations
[t,A] = ode45(f, 0:Dt:Time, IC); 

%% --- Data sampling ---
Time_head = 0;
Step_head = Time_head/Dt + 1;
Time_end = 5000;                   
Step_end = Time_end/Dt;
Time = Time_end - Time_head;        
data = A(Step_head:Step_end,:);     % samples of analysis 
window = 400;                       % time scale of analysis window 

%% --- Galerkin-based Granger causality (GGC) analysis ---
alpha = 0.1;        % significance level
permutation = 100;  % number of permutations for significance test (50-200),
                    % Np = 0: no permutation test
overlap = window/2;                     % time overlap between windows
num = floor((Time-window)/overlap+1);   % number of analysis windows
GGC_mat = zeros(num,9,9);               % GGC results
Pval_F = zeros(num,9,9);                % P value of F-test
Pval_LR = zeros(num,9,9);               % P value of  LR-test
GGC_perm = zeros(num,9,9);              % significance threshold calculated by permutation test
Corr_mat = zeros(num,9,9);              % Pearson corrlation
disp('computing GGC...')
tic();
% computing GGC in parallel
parfor i = 1:num
    % window segmentation
    Time_head = overlap*(i-1);
    Time_end = Time_head + window;
    Step_head = Time_head/Dt + 1;
    Step_end = Time_end/Dt;
    if Step_end > size(data,1)
       Step_end = size(data,1); 
    end
    GCGdata = data(Step_head:Step_end,:);
    % GGC analysis in windows are distributed across available CPU cores
    % If windows > CPU cores, computation / permutation proceeds in multiple rounds
    [GGC_mat(i,:,:),Pval_F(i,:,:),Pval_LR(i,:,:),GGC_perm(i,:,:)] ...
        = GalerkinGC(GCGdata,'permutations',permutation,'lag',1,'alpha',alpha); % GGC
    Corr_mat(i,:,:) = corr(GCGdata);
end
toc();

GGC = sum(GGC_mat,1)/num;       % average short-term causality to obtain long-term causality
GGC = reshape(GGC,[9 9]);
Pval_F = sum(Pval_F,1)/num;
Pval_F = reshape(Pval_F,[9 9]);
Sig_F = Pval_F < alpha;         % significant GGC determined by F test
Pval_LR = sum(Pval_LR,1)/num;
Pval_LR = reshape(Pval_LR,[9 9]);
Sig_LR = Pval_LR < alpha;       % significant GGC determined by LR test
GGC_perm = sum(GGC_perm,1)/num;
GGC_perm = reshape(GGC_perm,[9 9]);
Sig_perm = GGC > GGC_perm;      % significant GGC determined by permutation test
Corr = abs(sum(Corr_mat,1)/num);
Corr = reshape(Corr,[9 9]);     % Pearson correlation

CDI = std(GGC(find(GGC~=0)),[],'all') / mean(GGC(find(GGC~=0)),'all');  % Causal discriminatory index

% % plot GGC results
figure('Position', [100, 100, 1500, 500]); 
subplot(1,2,1); CausalMap(GGC)
title('GGC value','FontName','Times')
set(gca, 'Fontsize',24)
subplot(1,2,2); CausalMap(double(Sig_perm))
title(['Significant at \alpha=',num2str(alpha),', ',num2str(permutation),' permutations'],'FontName','Times')
set(gca, 'Fontsize',24)

% % plot F-test / LR-test result
% figure('Position', [100, 100, 1500, 500]); 
% subplot(1,3,1); CausalMap(abs(GGC))
% title('GGC value','FontName','Times')
% set(gca, 'Fontsize',24)
% subplot(1,3,2); CausalMap(double(Pval_LR))
% title('p-value (LR test)','FontName','Times')
% set(gca, 'Fontsize',24)
% subplot(1,3,3); CausalMap(double(Sig_F))
% title(['Significant at \alpha=',num2str(alpha)],'FontName','Times')
% set(gca, 'Fontsize',24)

%% --- Results Visualization ---
% % plot mode coefficients (in Matlab)
% figure('Position', [100, 100, 1200, 400]); 
% plot([Dt:Dt:Time],data);
% xlabel('physical time', 'FontSize',16,'FontWeight','bold');
% ylabel('mode coefficient', 'FontSize',16,'FontWeight','bold');
% title('Evolution of 9-mode shear flow', 'FontSize',20,'FontWeight','bold')
% grid on;
% legend('a1','a2','a3','a4','a5','a6','a7','a8','a9');

% plot flow contours (in Matlab)
figure; Shearflow3Dflow(data(1,:))

% % plot flow contours (in Tecplot)
% ShearflowTecplot(data(1,:),Dt,'shearflow_tecplot')

