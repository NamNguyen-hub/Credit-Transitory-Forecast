function oos_RMSE = oos_forecast_error_ar(y,p)
% Benjamin Wong
% Monash University
% July 2019
% Takes first 20 years to train the sample before calculating the one step
% ahead out of sample root mean square error for an AR(p)
%
%% INPUTS
%
%y                        Time series
%p                        lags
%
%% OUTPUTS
%oos_RMSE                   Out of sample forecast error

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

oos_resid = [];

%%
for jj = 80:T-1
    [AR_coeff] = olsvar([zeros(p,1);Y(1:jj,1)],p,1);
    Y_F = AR_coeff'*X(jj+1,:)';
    
    oos_resid = [oos_resid;Y(jj+1,:)-Y_F];
end

oos_RMSE = sqrt(mean(oos_resid.^2));

end