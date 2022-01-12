function oos_RMSE = BN_BVAR_oos_RMSE(y,p,lambda,target_variable)
% Benjamin Wong
% Monash University
% July 2019
% Companion .m file to BN_BVAR
% Takes first 20 years to train the sample before calculating the one step
% ahead out of sample root mean square error
%
%% INPUTS
%
%y                  Time series
%p                  number of lags
%lambda             Shrinkage hyper-parameter (see Banbura, Giannone and Reichlin, JAE, 2010)
%target_variable    the index (or column in y) of the target variable. Output growth in our application
%
%% OUTPUTS
%oos_RMSE                   The one step ahead pseudo out of sample forecast error

%% Preliminaries

[T N] = size(y);

%demean time series
y = y - repmat(mean(y),T,1);
%backcast data
y = [zeros(p,N);y];

K = p*N;        %number of parameters per equation

%Create Data Matrices
Y = y(p+1:end,:);     %Cut Away first p lags (the backcasted stuff)
X = [];

for i = 1:p
    Z = y(p+1-i:end-i,:);
    X = [X Z];
end

%% Set up matrices for variable with Signal-to-noise ratio prior


oos_resid = [];

for jj = 80:T-1
    
    
    %% Set up dummy observations
    
    small_sig2 =  ones(N,1);
    Stack_AR_coeff = zeros(p,N);
    
    for i = 1:N
        %     calculate variance for each equation to set prior
        %     by fitting an AR(4) per equation
        [~,small_sig2(i,1),~,~,~] = olsvar([zeros(p,1);Y(1:jj,i)],p,1);
    end
 
    Y_d = [zeros(N*p,N);
        diag(sqrt(small_sig2))];

    
    for i = 1:p
            Y_d((i-1)*N+1:i*N,:) = diag(Stack_AR_coeff(i,:)'.*repmat(i,N,1).*sqrt(small_sig2))/lambda;
    end
    
    X_d = [kron(diag(1:p),diag(sqrt(small_sig2)/lambda));
        zeros(N,K)];
   
    %% Do Least Squares to get posterior
    
    Y_star = [Y(1:jj,:);Y_d]; X_star = [X(1:jj,:); X_d];
    
    %Get VAR coefficients
    A_post = (X_star'*X_star)\(X_star'*Y_star);
    
    %One Step ahead Forecasts
    Y_F = X(jj+1,:)*A_post;
    
    oos_resid = [oos_resid;Y(jj+1,:)-Y_F];
    
end

oos_RMSE = sqrt(mean(oos_resid(:,target_variable).^2));

end