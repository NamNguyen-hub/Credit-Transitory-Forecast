function BN_cycle = multivariate_BN_ols(y,p)
% Benjamin Wong
% RBNZ
% August 2016
% Does Multivariate Beveridge Nelson Decomposition by demeaning, then
% backcasting. Estimates VAR or AR using least squares
% For more details, See Morley 2002, Economic Letters
%
%y              vector(s) of variables
%p              lags
%%
%estimate VAR
[T N] = size(y);

%demean
y = y - repmat(mean(y),T,1);
%backcast
y = [zeros(p,N);y];

%estimate
[A,~,~,~,X] = olsvar(y,p,1);

%Calculate BN cycle
X = [X(2:end,:);
    y(end,:) X(end,1:end-N)];

Companion = [A';
            eye(N*(p-1)) zeros(N*(p-1),N)];

bigeye = eye(N*p);

phi = -Companion*((bigeye-Companion)\eye(N*p));
BN_cycle = phi*X';
BN_cycle = BN_cycle(1:N,:)';

end

