%%  Galerkin-based Granger Causality (GGC) method
%   DESCRIPTION:
%     The Galerkin-based Granger Causality (GGC) is a modified Granger causality tailored to 
%     the structure of Galerkin models. It is specifically designed to investigate causal interactions 
%     in Galerkin-type and other data-driven reduced-order models of complex fluid flows.
%   USAGE:
%     [GGC, Pval_F, Pval_LR, GGC_perm] = GalerkinGC(data, 'lag', input1, 'alpha', input2, 'permutations', input3)
%   INPUTS:
%     data        - [T x N] matrix, T time points, N variables
%     'lag'       - (optional) time lag (default: 1)
%     'alpha'     - (optional) significance level (default: 0.05)
%     'permutations' - (optional) number of surrogate permutations (default: 0, for no permutation test)
%   OUTPUTS:
%     GGC         - [nvars x nvars] matrix, GGC(i,j) = Galerkin-based GC from mode j to mode i
%     Pval_F      - p-values from F-test
%     Pval_LR     - p-values from Likelihood Ratio test
%     GGC_perm    - significance thresholds of GGC value from permutation test
%   REFERENCE:
%     see §4 of [Wang, Y., Kou, J., Noack, B.R. et al. (2026). Causal analysis of a turbulent shear flow model]
%
function [GGC,Pval_F,Pval_LR,GGC_perm] = GalerkinGC(data,varargin) 
     [nobs,nvars] = size(data); 
     default_lag = 1;       % default: lag=1 should be used, unless time-lag of the target system is known
     default_alpha = 0.05;  % default: 0.05
     default_Np = 0;        % default: no permutation test
     p = inputParser;
     addRequired(p,'data');
     addParameter(p,'lag',default_lag);
     addParameter(p,'alpha',default_alpha);
     addParameter(p,'permutations',default_Np);
     parse(p,data,varargin{:});
     lag = p.Results.lag;
     alpha = p.Results.alpha;
     Np = p.Results.permutations;
     error = zeros(nvars,nvars,nobs-lag);
     reduced_model_indices = zeros(nvars,lag*nvars*(nvars-1)/2);
     indices_built = 0;
     Var = zeros(nvars);
     RSS = zeros(nvars);
     % REDUCED model : model of i excludes any term involving j
     % FULL model    : model of i includes all terms
     % i is the 'Granger result', j is the 'Granger cause'
     for i = 1:nvars                    % loop for 'Granger result'
        %% FULL Galerkin model
         X = zeros(lag*(nvars^2+3*nvars)/2+1,nobs-lag);     % Galerkin terms (regression terms)
         Y = zeros(1,nobs);                                 % 'Granger result'(regression objective)
         QuadTermVec = zeros(lag,nvars*(nvars+1)/2);
         QuadTermMat = zeros(lag,nvars,nvars);
         for k = 1:nobs-lag             % loop for time steps
             k2 = lag+k-1;
             for r = 1:lag              % loop for time lags
                 count = 1;
                 for p = 1:nvars        % loop for variable 1
                     for q = 1:nvars    % loop for variable 2
                         if q >= p
                            QuadTermVec(r,count) =  data(k+r-1,p) * data(k+r-1,q);
                            count = count + 1;
                            if (~indices_built) 
                                QuadTermMat(r,p,q) =  data(k+r-1,p) * data(k+r-1,q);
                                QuadTermMat(r,q,p) =  QuadTermMat(r,p,q);
                            end
                         end
                     end
                 end
             end
             X(:,k) = [reshape([data(k:k2, :),QuadTermVec],[lag*(nvars^2+3*nvars)/2,1]); 1]; % Galerkin terms
             % term indices for the REDUCED Galerkin model
             if (~indices_built)    % constructed only once (reused across iterations)
                 for p = 1:nvars
                     count = 1;
                     for q = 1:length(X(:,k))
                         if ismember(X(q,k), data(k:k2, p))             % remove linear terms involving 'cause' mode
                             continue;
                         elseif ~ismember(X(q,k), QuadTermMat(:,p,:))   % remove quadratic terms involving 'cause' mode
                             reduced_model_indices(p,count) = q;
                             count = count + 1;
                         end
                     end
                 end
             end
             if (~any(reduced_model_indices(:)==0))
                indices_built = 1;  % constructed 
             end
         end
         Y = data(:,i)';                            % time series of the 'Granger result'
         Af = Y(lag+1:end)/X;                       % regression of the full Galerkin model
         error(i,i,:) = Af*X - Y(lag+1:end);        % one-step fit, prediction error time series of the full model...
                                                    % for variable i are stored in error(i,i,:) 
         Var(i,i) = var(error(i,i,:));              % variance
         RSS(i,i) = sum(error(i,i,:).^2,'all');     % residual sum of squares (RSS)
        %% REDUCED Galerkin model
         for j = 1:nvars    % variable j is the 'Granger cause', which is removed in the reduced Galerkin model
             if i==j
                 continue
             end
             Y = zeros(1,nobs);
             Xr = X(reduced_model_indices(j,:),:);  % Galerkin terms excluding mode j
             Y = data(:,i)';    
             Ar = Y(lag+1:end)/Xr;                  % regression of the reduced model
             error(i,j,:) = Ar*Xr - Y(lag+1:end);   % one-step fit, prediction error time series removing variable j ...
                                                    % from the reduced model of variable i are stored in error(i,j,:) 
             Var(i,j) = var(error(i,j,:));              % variance
             RSS(i,j) = sum(error(i,j,:).^2,'all');     % residual sum of squares
         end
     end
     GGC = log(Var./diag(Var));         % Granger Causality: log-ratio of the variances from the two models
     %% Permutations test (recommanded for complex dynamics)
      Var_perm = zeros(Np,nvars,nvars);
      RSS_perm = zeros(Np,nvars,nvars);
      GGC_perm = zeros(nvars);
      shift = randi([100,nobs],1,Np);   % random time shift>100 of cause variable to break causality
      if Np == 0            % output 0 significance in case of no permutation test (infinity threshold)
          GGC_perm = inf * ( ones(nvars) - 2*diag(ones(nvars,1)) ); 
      end
      for n = 1:Np          % loop for permutations
         for j = 1:nvars    % variable j is the 'Granger cause', which is permutated
             data_perm = data;  
             data_perm(:, j) = circshift(data(:, j), shift(n)); % time shift using periodic boundary conditions
             X = zeros(lag*(nvars^2+3*nvars)/2+1,nobs-lag);
             QuadTermVec = zeros(lag,nvars*(nvars+1)/2);
             for k = 1:nobs-lag
                 k2 = lag+k-1;
                 for r = 1:lag
                     count = 1;
                     for p = 1:nvars
                         for q = 1:nvars
                             if q >= p
                                QuadTermVec(r,count) =  data_perm(k+r-1,p) * data_perm(k+r-1,q);
                                count = count + 1;
                             end
                         end
                     end
                 end
                 X(:,k) = [reshape([data_perm(k:k2, :),QuadTermVec],[lag*(nvars^2+3*nvars)/2,1]); 1];
             end
             for i = 1:nvars
                 if i == j
                     Var_perm(n,i,j) = Var(i,j);
                     RSS_perm(n,i,j) = RSS(i,j);
                     continue
                 end
                 Y = data_perm(:,i)';
                 Af_perm = Y(lag+1:end)/X;
                 error(i,j,:) = Af_perm*X - Y(lag+1:end);           % one-step fit, variable j is permutated
                 Var_perm(n,i,j) = var(error(i,j,:)); 
                 RSS_perm(n,i,j) = sum(error(i,j,:).^2,'all');
                 perms_GGC(n,i,j) = log(Var(i,j)/Var_perm(n,i,j));  % GC: log-ratio of the variances
             end
         end
         fprintf('Processing permutations for significance test... %d / %d\n', n, Np);
      end
      % result of permutation test
      if Np >= 1
         for i = 1:nvars
             for j = 1:nvars
                 if i == j
                    GGC_perm(i,j) = inf;    % output 0 significance on the diagonal (infinity threshold)
                    continue;
                 end
                GC_perm_ij = squeeze(perms_GGC(:,i,j));     % null distribution of hypothesis 'no GC from j to i'
                threshold = quantile(GC_perm_ij, 1-alpha);  % significance threshold is taken as (1-𝛼)-quantile of the distribution
                GGC_perm(i,j) = threshold;                  % output significance threshold
             end
         end
     end
     %% F-test (only applicable to simple dynamics)
     stat = zeros(nvars);
     df1 = size(X,1) - size(Xr,1);
     df2 = nobs - size(Xr,1) - 1;
     for i = 1:nvars
         for j = 1:nvars
             stat(i,j) = df2 / df1 * (RSS(i,j)/RSS(i,i)-1); 
         end
     end
     Pval_F = 1 - fcdf(stat, df1, df2);         % output p value
     %% Likelihood-ratio test (only applicable to simple dynamics)
     stat = zeros(nvars);
     df1 = size(X,1) - size(Xr,1);
     for i = 1:nvars
         for j = 1:nvars
             stat(i,j) = log(RSS(i,j)/RSS(i,i)); 
         end
     end
     Pval_LR = 1 - chi2cdf((nobs-1)*stat,df1);  % output p value
end

