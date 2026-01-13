%%  Galerkin-based GC analysis on turbulent lid-driven cavity flow
%   DISCRIPTION:
%     Perform causal analysis on the dominant POD modes of a turbulent cavity flow, utilizing proposed GGC method.
%     See appendix C / supplymentary material of
%     [Wang, Y., Kou, J., Noack, B.R. et al. (2026). Causal analysis of a turbulent shear flow model]
%
%   DATA SOURCE:
%     [Arbabi, H., & Mezić, I. (2017). Study of dynamics in post-transient flows using Koopman mode decomposition]
%     flow data can be found at: https://mgroup.me.ucsb.edu/resources 
%
clear all
addpath('./Data_CavityFlow');
addpath('./Code_CavityFlow');
addpath('./Code_SINDy');
addpath('./Code_PlotResult');

%% --- Load cavity flow data ---
num_snap = 10000;

% % --- Approach 1. load flow data and perform POD (slow) ---
% % Flow data is not provided, which can be found at https://mgroup.me.ucsb.edu/resources.
% % This process can be slow, an alternative is to directly load POD modes through Approach 2.
% disp('Loading flow data...')
% [VelField] = compute_velfield('Cavity30k.mat'); 
% Dt = VelField.t;
% Dt = Dt(2);
% Nodes = length(VelField.xv)^2;
% UX = VelField.UX;
% UX = UX(:,:,1:num_snap);
% UX = reshape(UX, Nodes, num_snap);
% UY = VelField.UY;
% UY = UY(:,:,1:num_snap);
% UY = reshape(UY, Nodes, num_snap);
% gridx = VelField.xv;                      
% gridy = VelField.yv;    
% % perform POD on velocity
% Snap = [UX; UY];    % concatenating the u- and v-velocity components into a unified matrix.
% Snap = bsxfun(@minus,Snap,mean(Snap,2));  % demean
% disp('Performing POD...')
% [U,S,V] = svd(Snap(:,1:end),'econ');      % POD decomposition 
% mode_pod = U;                         	% POD modes
% coef_pod = U' * Snap;                    	% POD coefficients
% energy_pod = diag(S);                    	% mode energy             

% --- Approach 2. load POD data directly (fast) ---
load Cavity16k_POD50.mat    % only first 50 POD modes
disp('POD modes loaded')
mode_pod = pod.mode;
coef_pod = pod.coef;
energy_pod = pod.energy;
Dt = pod.dt;
gridx = pod.gridx;
gridy = pod.gridy;
Nodes = (length(gridx))^2;

% % --- plot energy spectrum of POD modes ---
% norm_energy = energy_pod / sum(energy_pod);
% PODEnergyPlot(norm_energy);       

%% --- Galerkin-based Granger causality (GGC) analysis ---
num_mode = 10;      % truncation of POD modes
mode_pod = mode_pod(:,1:num_mode);
data = coef_pod(1:num_mode,:)';
Time_head = 0;
Step_head = Time_head/Dt + 1;
Time_end = num_snap*Dt;
Step_end = num_snap;
Time = Time_end - Time_head;
Steps = Step_end - Step_head + 1;

window = 5000*Dt;                           % length of analysis window (in physical time)
overlap = window/2;                         % time overlap between windows
num = floor((Time-window)/overlap+1);       % number of analysis windows

alpha = 0.05;       % significance level
test = 'F';         % type of statstic test: 'F' for F-test, 'LR' for log Likelihood-ratio test
                    % F-test and LR test are efficient, but much less reliable compared to permutation test
permutation = 100;  % number of permutations for significance test
                    % permutation = 0: no permutation test
GGC_mat = zeros(num,num_mode,num_mode);     % GGC results
Pval_F = zeros(num,num_mode,num_mode);      % P values of F-test
Pval_LR = zeros(num,num_mode,num_mode);     % P values of LR-test
GGC_perm = zeros(num,num_mode,num_mode);    % significance threshold calculated by permutation test
Corr_mat = zeros(num,num_mode,num_mode);    % Pearson corrlation
% computing GGC in parallel
disp('Computing GGC...')
tic();
parfor i = 1:num
    % window segmentation
    Time_head = overlap*(i-1);
    Time_end = Time_head + window;
    Step_head = Time_head/Dt + 1;
    Step_end = Time_end/Dt;
    if Step_end > size(data,1)
       Step_end = size(data,1); 
    end
    GCdata = data(Step_head:Step_end,:);
    % GGC analysis in windows are distributed across available CPU cores
    % If windows > CPU cores, computation / permutation proceeds in multiple rounds
    [GGC_mat(i,:,:),Pval_F(i,:,:),Pval_LR(i,:,:),GGC_perm(i,:,:)] ...
        = GalerkinGC(GCdata,'permutations',permutation,'lag',1,'alpha',alpha);  % GGC
end
toc();

GGC = sum(GGC_mat,1)/num;       % average short-term causality to obtain long-term causality
GGC = reshape(GGC,[num_mode num_mode]);
Pval_F = sum(Pval_F,1)/num;
Pval_F = reshape(Pval_F,[num_mode num_mode]);
Sig_F = Pval_F < alpha;         % significant GGC determined by F test
Pval_LR = sum(Pval_LR,1)/num;
Pval_LR = reshape(Pval_LR,[num_mode num_mode]);
Sig_LR = Pval_LR < alpha;       % significant GGC determined by LR test
GGC_perm = sum(GGC_perm,1)/num;
GGC_perm = reshape(GGC_perm,[num_mode num_mode]);
Sig_perm = GGC > GGC_perm;      % significant GGC determined by permutation test

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
% subplot(1,3,2); CausalMap(double(Pval_F))
% title('p-value (F test)','FontName','Times')
% set(gca, 'Fontsize',24)
% subplot(1,3,3); CausalMap(double(Sig_F))
% title(['Significant at \alpha=',num2str(alpha)],'FontName','Times')
% set(gca, 'Fontsize',24)

%% --- SINDy regression for Sparse Galerkin model (optional) ---
% % Optional: Sparse regression of Galerkin model via SINDy for validation (see supplement material)      
% n = 4;
% dx = zeros(num_snap, n);
% for i = 1:n
%     dx(:, i) = gradient(coef_pod(i, :), Dt)';
% end
% nosine = 0;
% polyorder = 2;
% Theta = poolData(coef_pod(1:n, :)', n, polyorder, nosine);
% lambda = 0.1;
% Xi = sparsifyDynamics(Theta, dx, lambda, polyorder);
% cell_list = arrayfun(@(x) sprintf('X%d', x), 1:n, 'UniformOutput', false);
% poolDataLIST(cell_list, Xi, n, polyorder, nosine);
% t_end = 30;
% options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8 * ones(n, 1));
% tspan = [Dt:Dt:t_end];
% x0 = coef_pod(1:n, 1);  
% [tt, xt] = ode45(@(t,x) sparseGalerkin(t, x, Xi, polyorder, nosine), tspan, x0, options);
% tend_index = round(t_end / Dt);  
% % Create combined comparison plot 
% figure('Position', [200, 200, 1400, 500]);
% colors = lines(n);
% line_styles_original = {'-', '-', '-', '-'};
% line_styles_fitted = {'--', '--', '--', '--'};
% hold on;
% for i = 1:n
%     plot(tt, coef_pod(i, 1:tend_index)', 'Color', colors(i, :), 'LineStyle', line_styles_original{i}, ...
%          'LineWidth', 3.0, 'DisplayName', sprintf('X%d Original', i));
%     plot(tt, xt(:, i), 'Color', colors(i, :), 'LineStyle', line_styles_fitted{i}, ...
%          'LineWidth', 2.5, 'DisplayName', sprintf('X%d Fitted', i));
% end
% xlabel('Time', 'FontWeight', 'bold');
% ylabel('Mode Amplitude', 'FontWeight', 'bold');
% legend('Location', 'best', 'NumColumns', 2);
% set(gca, 'GridAlpha', 0.4, 'MinorGridAlpha', 0.2);
% set(gca, 'FontSize', 20); 
% set(gca, 'FontName','times'); 
% set(gcf, 'Color', 'white');
% xlim([0 35])
% hold off;

%% --- Results Visualization ---
% plot POD modes (in Matlab)
Vel = sqrt(mode_pod(1:Nodes,:).^2 + mode_pod(Nodes+1:end,:).^2);    % velocity magnitude
PlotVorticity(Vel(:,1),1);

% % plot POD modes (in Tecplot format)
% [x,y] = meshgrid(gridx,gridy);                        
% x = reshape(x,[],1);            % grid x      
% y = reshape(y,[],1);        	% grid y   
% Vel = sqrt(mode_pod(1:Nodes,:).^2 + mode_pod(Nodes+1:end,:).^2);    % velocity magnitude
% DataCell = {'X',x, 'Y',y };
% for i = 1:num_mode
%     DataCell = [ DataCell, ['Ux',num2str(i)] ,mode_pod(1:Nodes,i), ...
%                 ['Uy',num2str(i)], mode_pod(Nodes+1:end,i), ...
%                 ['Vel',num2str(i)], Vel(:,i) ];
% end
% CavityTecplot('Cavity16k_POD', DataCell, sqrt(Nodes), sqrt(Nodes));

% % plot POD coefficients (in Matlab) 
% figure('Position', [100, 100, 1200, 400]); 
% plot(coef_pod'); xlim([1 1000]); 
% xlabel('time','FontWeight','bold');
% ylabel('mode coefficient','FontWeight','bold');
% set(gca,'FontSize',16)

% % plot POD coefficients (in Tecplot format)
% TecplotLines(coef_pod(1:num_mode,:)', Dt, 'Cavity16k_PODcoef.dat')




