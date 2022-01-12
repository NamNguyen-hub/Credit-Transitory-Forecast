%% Pseudo_real_time

PRT_gap = cell(size(dates,1),1);

parpool(4)

lambda_store = NaN(size(dates,1),1);

lambda0 = 1;
lambdafull = fminsearch(@(lambda0) BN_BVAR_oos_RMSE(y{2},p,lambda0,target_variable(2)),0,options.optimisation);
[BN_cycle_full] = BN_BVAR(y{2},p,lambdafull,target_variable(2));
Final_gap = [dates BN_cycle_full(:,target_variable(2))];

parfor jj = 120:size(dates,1)
    
    %% Minimising One step ahead OOS RMSE
    warning('off', 'MATLAB:illConditionedMatrix');
    lambda0 = 1;
    lambda = fminsearch(@(lambda0) BN_BVAR_oos_RMSE(y{2}(jj-120+1:jj,:),p,lambda0,target_variable(2)),0,options.optimisation);
    
    
    %%
    [BN_cycle] = BN_BVAR(y{2}(jj-120+1:jj,:),p,lambda,target_variable(2));
    [BN_cycle_lambdafull] = BN_BVAR(y{2}(jj-120+1:jj,:),p,lambdafull,target_variable(2));
    
    temp_BN = [BN_cycle(:,target_variable(2)) BN_cycle_lambdafull(:,target_variable(2))];
    RT_gap_temp = [dates(jj-120+1:jj) temp_BN];
    
    PRT_gap{jj,1} = RT_gap_temp;
    lambda_store(jj,1) = lambda;
end

PRT_gap(1:119,:) = [];

delete(gcp)

%%

for jj = 1:size(PRT_gap,1)
    
    RT_gap(jj,:) = PRT_gap{jj,1}(end,:);
    
end

%%
figure
subplot(2,1,1)
h1=NBERbc(Final_gap(120:end,1),Final_gap(120:end,2),{'-'},3,{'r'});
hold on
h2=plot(RT_gap(:,1),RT_gap(:,2),'-.b','LineWidth',3);
hold on
plot([RT_gap(1,1) RT_gap(end,1)],zeros(2,1),'-k','LineWidth',1);
set(gca,'FontSize',14)
title('Reoptimizing \lambda in Pseudo Real Time')
ylim([-6 4])
legend([h1 h2],{'Benchmark','Rolling'},'Orientation','horizontal')

subplot(2,1,2)
h1=NBERbc(Final_gap(120:end,1),Final_gap(120:end,2),{'-'},3,{'r'});
hold on
h2=plot(RT_gap(:,1),RT_gap(:,3),'-.b','LineWidth',3);
hold on
plot([RT_gap(1,1) RT_gap(end,1)],zeros(2,1),'-k','LineWidth',1);
set(gca,'FontSize',14)
title('Keeping \lambda Fixed in Pseudo Real Time')
ylim([-6 4])
legend([h1 h2],{'Benchmark','Rolling'},'Orientation','horizontal')

disp('Figure 9 Done')
toc
