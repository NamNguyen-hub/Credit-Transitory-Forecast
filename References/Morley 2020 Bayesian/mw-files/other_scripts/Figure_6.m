%% DIFFERENT INFORMATION SETS
breakpoint = find(dates ==2006.25);
for jj = 1:3
    
    % Demeaning relative to break
    
    T1 = size(y{jj}(1:breakpoint,:),1);
    T2 = size(y{jj}(breakpoint+1:end,:),1);
    
    mean1 = mean(y{jj}(1:breakpoint,:));
    mean2 = mean(y{jj}(breakpoint+1:end,:));
    
    y_b = [y{jj}(1:breakpoint,:)-repmat(mean1,T1,1);...
        y{jj}(breakpoint+1:end,:)-repmat(mean2,T2,1)];
    
    %%

    lambda_break = fminsearch(@(lambda0) BN_BVAR_oos_RMSE(y_b,p,lambda0,target_variable(jj)),0,options.optimisation);
    [BN_cycle_break] = BN_BVAR(y_b,p,lambda_break,target_variable(jj));
    
    Output_gap_breaks(:,jj) = BN_cycle_break(:,target_variable(jj));
    
    
end

%%
% Run Benchmark here
lambda = fminsearch(@(lambda0) BN_BVAR_oos_RMSE(y{2},p,lambda0,target_variable(2)),0,options.optimisation);
[BN_cycle] = BN_BVAR(y{2},p,lambda,target_variable(2));

%% Plots
figure
subplot(2,1,1)
h1=NBERbc(dates,BN_cycle(:,target_variable(2)),{'-'},3,{'r'});
hold on
h2=plot(dates,Output_gap_breaks(:,2),'-.k','LineWidth',3);
hold on
plot([dates(1) dates(end)],zeros(2,1),'-k','LineWidth',2)
legend([h1 h2],{'Benchmark','Allowing for Break in Drift'},'Location','southeast','Orientation','Horizontal')
set(gca,'FontSize',16)
title('Estimated Output Gap from 23 Variable Benchmark')


subplot(2,1,2)
h1=NBERbc(dates,Output_gap_breaks(:,1),{':'},3,{'b'});
hold on
h2=plot(dates,Output_gap_breaks(:,2),'-.r','LineWidth',3);
hold on
h3=plot(dates,Output_gap_breaks(:,3),'-','LineWidth',3);
hold on
plot([dates(1) dates(end)],zeros(2,1),'-k','LineWidth',2);
legend([h1 h2 h3],{'8 Variables','23 Variables','138 Variables'},'Location','southeast','Orientation','Horizontal')
set(gca,'FontSize',16)
title('Allowing for Break in Drift in 2006Q2')


disp('Figure 6 Done')
toc

