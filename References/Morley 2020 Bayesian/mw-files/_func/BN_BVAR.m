function [BN_cycle,Info_decom,shock_decom,FEVD,Std_err] =...
    BN_BVAR(y,p,lambda,target_variable,varargin)
% Benjamin Wong
% Monash University
% July 2019
% Estimates large BVAR by using the Normal-Wishart version of the Minnesota Prior.
% Subsequently performs BN decomposition and decomposes the BN trend and cycle into sources of multivariate
% information and shocks like Morley and Wong (JAE)
% Identification by Cholesky Decomposition
%
%
%% INPUTS
%
%y                  Time series
%p                  number of lags
%lambda             Shrinkage hyper-parameter (see Banbura, Giannone and Reichlin, JAE, 2010)
%target_variable    the index (or column in y) of the target variable. Output growth in our application
%varargin           input 'Decomposition' to do information and structural decomposition
%                   input ' Standard_Errors' to generate BN cycle standard errors.
%                   This calculation is documented in Section A4 of the
%                   online appendix of Kamber, Morley and Wong (REStat,
%                   2018)
%
%% OUTPUTS
%BN_cycle                   The estimated BN cycle
%
% The following outputs will only be produced with the varargin options
%
%Info_decom                 Stores the informational decomposition of both trend and
%                           cycle
%shock_decom                Stores the structural decomposition of both
%                           trend and cycle
%FEVD                       Forecast error variance decomposition of trend
%                           and cycle.
%Std_err                    Reports the standard errors of the BN cycle
%                           according to the calculation in Kamber, Morley
%                           and Wong (REStat, 2018)

%% Preliminaries

[T N] = size(y);

%demean time series
y = y - repmat(mean(y),T,1);
%backcast data
y = [zeros(p,N);y];

K = p*N;        %number of parameters per equation

small_sig2 =  zeros(N,1);
Stack_AR_coeff = zeros(p,N);

for i = 1:N
    %     calculate variance for each equation to set prior
    %     by fitting an AR(4) per equation
    [~,small_sig2(i,1),~,~,~] = olsvar(y(:,i),p,1);
end

%Create Data Matrices
Y = y(p+1:end,:);     %Cut Away first p lags (the backcasted stuff)
X = [];

for i = 1:p
    Z = y(p+1-i:end-i,:);
    X = [X Z];
end

%% Set up dummy observations


Y_d = [zeros(N*p,N);
    diag(sqrt(small_sig2))];

for i = 1:p
    Y_d((i-1)*N+1:i*N,:) = diag(Stack_AR_coeff(i,:)'.*repmat(i,N,1).*sqrt(small_sig2))/lambda;
end

X_d = [kron(diag(1:p),diag(sqrt(small_sig2)/lambda));
    zeros(N,K)];

%% Do Least Squares to get posterior
Y_star = [Y;Y_d]; X_star = [X; X_d];

%Get VAR coefficients
A_post = (X_star'*X_star)\(X_star'*Y_star);

%Get BVAR residuals and calcuate posterior covariance matrix
U = Y-X*A_post;
U_star = Y_star-X_star*A_post;
SIGMA = (U_star'*U_star)/(size(Y_star,1)-size(A_post,1));


%% BN Decomposition Starts here
%Add current period in to incoporate current information set to do BN decomposition

X = [X(2:end,:);
    Y(end,:) X(end,1:end-N)];

Companion = [A_post';
    eye(N*(p-1)) zeros(N*(p-1),N)];

bigeye = eye(N*p);

phi = -Companion*((bigeye-Companion)\eye(N*p));
BN_cycle = phi*X';
BN_cycle = BN_cycle(1:N,:)';

%% Do decompositons of BN trend and cycle
% Setup inputs to calculating the informational and shock decomposition
bigeye_sec = [eye(N) zeros(N,N*(p-1))];
Info_decom = [];
shock_decom = [];
FEVD = [];
Std_err=[];
selector_vec = zeros(1,N);
selector_vec(target_variable) = 1;


if sum(strcmp(varargin,'Decomposition')) == 1
    %construct historical structural shocks
    
    eta = zeros(size(U));
    A0 = chol(SIGMA,'lower');
    A0_inv = A0\eye(N);
    
    for ii = 1:T
        eta(ii,:) = A0_inv*U(ii,:)';
        
    end
    
    Info_decom.cycle = zeros(size(X,1),N);
    Info_decom.trend = zeros(size(X,1),N);
    
    shock_decom.cycle = zeros(size(X,1),N);
    shock_decom.trend = zeros(size(X,1),N);
    
    % see formulaes in the paper
    for jj = 1:size(X,1)
        count = 0;
        
        for ii = jj:-1:1
            shock_vec = diag(U(ii,:));
            s_shocks_vec = diag(eta(ii,:));
            if ii == jj
                Info_decom.trend(jj,:) = selector_vec*bigeye_sec*((bigeye-Companion)\eye(N*p))*bigeye_sec'*shock_vec;
                shock_decom.trend(jj,:) = (selector_vec*bigeye_sec*((bigeye-Companion)\eye(N*p))*bigeye_sec'*A0*s_shocks_vec)';
            end
            
            Info_decom.cycle(jj,:)= Info_decom.cycle(jj,:)+ ...
                selector_vec*bigeye_sec*(phi*(Companion^count))*bigeye_sec'*shock_vec;
            
            shock_decom.cycle(jj,:)= shock_decom.cycle(jj,:)+ ...
                selector_vec*bigeye_sec*(phi*(Companion^count))*bigeye_sec'*A0*s_shocks_vec;
            
            count = count + 1;
        end
    end
    %% Calculate FEVD
    
    nFEVD = 101;    %number of FEVD periods to compute
    FEVD.cycle = zeros(N,nFEVD); %Group by variables
    
    for i = 1:nFEVD
        temp = bigeye_sec*phi*(Companion^(i-1))*bigeye_sec'*A0;
        temp = temp(target_variable,:).^2;
        if i == 1
            FEVD.cycle(:,1) = temp;
        else
            FEVD.cycle(:,i) = FEVD.cycle(:,i-1)+temp';
        end
    end
    %trend
    BigPhi = (bigeye-Companion)\eye(N*p);
    FEVD.trend = reshape(bigeye_sec*BigPhi*bigeye_sec'*A0,[],N);
    FEVD.trend = (FEVD.trend(target_variable,:).^2)';
    FEVD.trend = 100*FEVD.trend./repmat(sum(FEVD.trend),N,1);
    
end

%% Calculate Standard  Error

if sum(strcmp(varargin,'Standard_Errors')) == 1
    Q = zeros(size(Companion));
    Q(1:N,1:N) = SIGMA;
    
    % See Kamber Morley and Wong (REStat,2018, online appendx)
    SIGMA_X = (eye(size(Companion,1)^2)-kron(Companion,Companion))\Q(:);
    SIGMA_X = reshape(SIGMA_X,N*p,N*p);
    Std_err = phi*SIGMA_X*phi';
    Std_err = sqrt(diag(Std_err(1:N,1:N)));
    
end



end