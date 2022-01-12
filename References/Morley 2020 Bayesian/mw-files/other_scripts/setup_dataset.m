% We just turn off this warning because it can a bit annoying. Especially
% within the optimization portion of the procedure.

% Set up other things here
warning('off', 'MATLAB:illConditionedMatrix');

% Optimization options and starting value for lambda
options.optimisation = optimset('Display','off','TolX',1e-8);
lambda0 = 1;

%================================================================================
%% SET UP DATASET
%================================================================================
% There are three code for the VAR_Size
%VAR_SIZE = 2;     %1-BASIC, 2-BASELINE, 3-LARGE

%Position of out output growth in the 3 datasets. Output growth is the target variable as we are estimating output gap
target_variable = [1;2;2];

rawdata.FRED = xlsread('multivariate_BN_dataset.xlsx','FREDQD','B7:IT286');
rawdata.haver_series = xlsread('multivariate_BN_dataset.xlsx','haver_series','B4:C283');

% Code 1 is slow moving 2 is fast moving. This distinction is only for the
% structural analysis
processing_code = xlsread('multivariate_BN_dataset.xlsx','FRED_Mnenomic','B2:G258');

% These just list the series we use for the benchmark model
Series = {'Oil Price','GDP growth','PCE','Income','IP','CAPU','Employment','Unemployment','Hours','Housing Starts','PCE inflation','GDP Deflator','CPI','PPI','Hourly Earnings','Productivity','Fed Funds Rate','Yield Curve','Real M1','Real M2','TR','NBR','Stock Price'};

for jj = 1:3
    VAR_SIZE = jj;
    
    % This distinction is just for the structural identification ordering
    slow_fast_moving_code = [find(processing_code(:,VAR_SIZE+1)==1);find(processing_code(:,VAR_SIZE+1)==2)];
    y{jj,1} = rawdata.FRED;
    
    % Just for nonborrowed reserves because it turns negative. We do
    % percentage change instead of log first difference
    y{jj,1}(:,209) = [NaN;100*(y{jj,1}(2:end,209)./y{jj,1}(1:end-1,209)-1)];
    
    % Pulling out datasets
    y{jj,1} = y{jj,1}(:,slow_fast_moving_code);
    
    
    %Processing levels as logs
    
    log_code = zeros(size(y{jj,1},2),1); %Whether to log or not
    level_code=find(processing_code(slow_fast_moving_code,1));
    y{jj,1}(:,level_code) = 100*log(y{jj,1}(:,level_code));
    log_code(level_code,1) = 1;
    
    % Put oil price first in 23 and 138 variables. Drop oil in 8 variable version
    if VAR_SIZE > 1
        y{jj,1} =[100*log(rawdata.haver_series(:,2)) y{jj,1} 100*log(rawdata.haver_series(:,1))]; %add in oil and stock prices
        log_code = [1;log_code;1];
    else
        y{jj,1} =[y{jj,1} 100*log(rawdata.haver_series(:,1))]; %add in stock prices
        log_code = [log_code;1];
    end
    
    y{jj,1} = y{jj,1}(50:end,:); %Starts in 1959Q2
    
    % This is just for the break test for difference in mean for the first and
    % second half
    T = size(y{jj,1},1);
    t = round(T/2);
    
    % Get rid of series without first observation in 1959Q2 and with NaNs
    nan_variables = find(sum(isnan(y{jj,1})));
    y{jj,1}(:,nan_variables) = [];
    log_code(nan_variables,:) = [];
    
    N = size(y{jj,1},2);
    %output_level = y{jj,1}(2:end,2);
    
    Trans_code = zeros(N,1);
    % Do STATIONARITY ADJUSTMENT with t test
    for i = 1:N
        [h pvalue] = ttest2(y{jj,1}(1:t,i),y{jj,1}(t+1:end,i));
        if pvalue < 0.05
            y{jj,1}(:,i) = [NaN;diff(y{jj,1}(:,i))];
            Trans_code(i,1) = Trans_code(i,1)+1;
        end
        
    end
    
    y{jj,1} = y{jj,1}(2:end,:);
    
    % Do unit root test
    for i = 1:size(y{jj,1},2)
        adf= unitroot(y{jj,1}(:,i));
        if adf(3,4)>0.05
            y{jj,1}(:,i) = [0;diff(y{jj,1}(:,i))];
            Trans_code(i,1) = Trans_code(i,1)+1;
        end
    end
end