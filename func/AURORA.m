





function  [Outputs,running_time]= AURORA(X,opts)

% This implementation is to detect anomaly points in multivariate time
% series
% Linzhang @UAlbany
% collabrate with Wenyu Zhang @Cornell

% for period dictionary
if ~isfield(opts, 'Dictionary_type'),   opts.Dictionary_type   = 'Ramanujan'; end  %
if ~isfield(opts, 'Pmax'),                  opts.Pmax   = 50; end  %
if ~isfield(opts, 'R'),                 opts.R      = Create_Dictionary(opts.Pmax,size(X,1),opts.Dictionary_type); end  %
if ~isfield(opts, 'Penalty_type'),      opts.Penalty_type   = 'square'; end

% for B-spline dictionaty
if ~isfield(opts, 'D'),                 opts.D = difference_operator_matrix(opts.A(:,2:end),'third'); end  %
if ~isfield(opts, 'halting_Thr'),  opts.halting_Thr = 1e-3; end  % for algo
if ~isfield(opts, 'max_iter'),     opts.max_iter = 100; end  % for algo

if ~isfield(opts, 'visual'),       opts.visual = 0; end  % for algo

if ~isfield(opts, 'rho'),          opts.rho =1.1; end
if ~isfield(opts, 'is_show'),      opts.is_show =0; end


if ~isfield(opts, 'lambda_1'),      opts.lambda_1 = 0.1; end
if ~isfield(opts, 'lambda_2'),      opts.lambda_2 = 0.1; end
if ~isfield(opts, 'lambda_3'),      opts.lambda_3 = 0.1; end


Pmax        = opts.Pmax;
halting_Thr = opts.halting_Thr; % for outer loop
max_iter    = opts.max_iter;

lambda_1 = opts.lambda_1;
lambda_2 = opts.lambda_2;
lambda_3 = opts.lambda_3;

%% get period penalty
[periods,penalty_vector]  =  get_period_penalty(Pmax,opts.Penalty_type);
H_inv                     =  diag(1./penalty_vector); % inverse of diagnal matrix


%% initialization

G       = opts.R* H_inv;  % H is diagnal penlty matrix
A       = opts.A;         % B-spline dictionary

% intercept in spline as zero in D
D1         = opts.D; % difference operator matrix
[a,b]      = size(D1);
D          = zeros(a+1,b+1);
D(2:a+1,2:b+1) = D1;


rho_1  = 1;
rho_2  = 1;
rho_3  = 1;


[T,N]  = size(X);
[~,l1] = size(G);
[~,l2] = size(A);

U =  zeros(l1,N);
W =  zeros(l2,N);

Gamma_1 = zeros(size(U));
Gamma_2 = zeros(size(W));
Gamma_3 = zeros(size(X));


I_U =  eye(size(G,2));
I_W =  eye(size(A,2));

%% PRE-compute
GTX =  G'*X;
GTA =  G'*A;
DTD =  D'*D;
ATA =  A'*A;
ATG =  A'*G;
GTG =  G'*G;


QUIET    = 0;
t_start  = tic;
%% main loop


for iter = 1: max_iter
    
    
    % Update Y
    E = X- G*U -A*W -Gamma_3/rho_3;
    Y = sign(E).*max(abs(E) - 1/rho_3,0);
    
    % Update V
    H = U - Gamma_1/rho_1;
    V = sign(H).*max(abs(H)-lambda_1/rho_1,0);
    
    
    % Update P
    M           = W - Gamma_2/rho_2;
    [J,Sigma,K] = svd(M,'econ');
    Sigma_thr   = diag(max(diag(Sigma) - lambda_2/rho_2 ,0));
    P           = J*Sigma_thr*K';
    
    % Update U
    Tmp_U_1 = rho_1 * I_U + rho_3* GTG;
    Tmp_U_2 = rho_1* V + Gamma_1 + rho_3*(GTX - GTA*W - G'*Y) - G'*Gamma_3;
    U       = inv(Tmp_U_1) *Tmp_U_2;
    
    % Update W
    Tmp_W_1 = lambda_3*DTD  + rho_1 * I_W + rho_3* ATA;
    Tmp_W_2 = rho_2*P + Gamma_2 + rho_3* ATG*U - rho_3* A'*Y- A'*Gamma_3;
    W       = inv(Tmp_W_1) * Tmp_W_2;
    
    
    % Update Lagrangian multipliers
    Gamma_1  = Gamma_1 + rho_1 *(V - U);
    Gamma_2  = Gamma_2 + rho_2 *(P - W);
    Gamma_3  = Gamma_3 + rho_3 *(Y - (X- G*U -A*W));
    
    
    % ++++++++++++++++++++++ update penalty parameters
    %     rho_1  =  min(tau*rho_1,1e10); % using this min is bad.
    %     rho_2  =  min(tau*rho_2,1e10);
    %     rho_3  =  min(tau*rho_3,1e10);
    
    %++++++++++++++++++++++ check convergence
    if iter > 1
        history.objval(iter)  = getObj(X,G,U,A,W,D,opts);
        
        obj_resi =  history.objval(iter)- history.objval(iter-1);
        % obj_resi              = abs(history.objval(iter)- history.objval(iter-1))/history.objval(iter-1);
        if mod(iter,5) == 0
            disp(['Iter:',num2str(iter),'::','obj_res = ',num2str(abs(obj_resi))]);
        end
        
        if obj_resi <  halting_Thr
            break;
        end
    end
    
end




%% output

All_variables.Y = Y;
All_variables.V = V;
All_variables.P = P;
All_variables.U = U;
All_variables.W = W;
All_variables.G = G;



Outputs.All_variables = All_variables;

% if ~QUIET
%     running_time = toc(t_start);
% end

running_time = toc(t_start);
%% GET PERIOD
% detect_periods_vector = cell(1,N);
for kk =  1:N
    s_individual          =  H_inv * V(:,kk);
    periods_vector_tmp    =  compute_vector_period(s_individual,opts);
    
    [periods_output_mag, periods_idx_tmp] = sort(periods_vector_tmp,'descend');
    detected_periods(kk,:)                = periods_idx_tmp;% (1:TOP_K); % top-1
    
end

Outputs.detected_periods = detected_periods;



%% GET outlier

detected_trends_data_matrix   =  A*P;
detected_periodic_data_matrix =  G*V;

Outlier = X - detected_trends_data_matrix - detected_periodic_data_matrix;
%Outlier =  All_variables.Z;

for jj =  1:size(Outlier,2)
    [outliers_output_mag, outliers_idx_tmp] = sort(Outlier(:,jj),'descend');
    detected_outliers(jj,:)                 = outliers_idx_tmp;% (1:TOP_K); % top-1
end

Outputs.detected_outliers = detected_outliers;

Outputs.O = Outlier;
end




function obj = getObj(X,G,U,A,W,D,opts)
obj =  norm(X-G*U-A*W,1) + opts.lambda_1*norm(W,1) +opts.lambda_2*nuclear_norm(U) +opts.lambda_3*norm(D*W,'fro');
end




function detect_periods_vector = compute_vector_period(s,opts)

Pmax              = opts.Pmax;
energy_s          = 0.*[1:Pmax];
current_index_end = 0;

for i=1:Pmax
    k_orig = 1:i;k=k_orig(gcd(k_orig,i)==1);
    current_index_start = current_index_end + 1;
    current_index_end   = current_index_end + size(k,2);
    
    for j = current_index_start:current_index_end
        energy_s(i) = energy_s(i)+((abs(s(j)))^2);
    end
end

energy_s(1)    = 0;
detect_periods_vector =  energy_s./norm(energy_s);
end
