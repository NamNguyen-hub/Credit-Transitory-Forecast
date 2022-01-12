% RUN Benchmark again
lambda = fminsearch(@(lambda0) BN_BVAR_oos_RMSE(y{2},p,lambda0,target_variable(2)),0,options.optimisation);
[BN_cycle,Info_decom,shock_decom,FEVD] = BN_BVAR(y{2},p,lambda,target_variable(2),'Decomposition');

Totals = [sum(sum(Info_decom.trend(83:95,:))) sum(sum(Info_decom.trend(195:218,:)))];

%%

figure
subplot(2,1,1)
bar([sum(Info_decom.trend(83:95,3)) sum(Info_decom.trend(195:218,3))],'k');
xlim([0.5 2.5])
set(gca,'xticklabel',{'1980Q1-1983Q1','2008Q1-2013Q4'})
set(gca,'FontSize',16)
title('Cumulative Change in Trend Accounted for by Consumption Forecast Errors')
subplot(2,2,3)
NBERbc(dates,y{2}(:,3),{'-'},3,{'r'});
hold on
plot([dates(1) dates(end)],zeros(2,1))
xlim([1980 1983])
set(gca,'FontSize',16)
title('Real PCE')
subplot(2,2,4)
NBERbc(dates,y{2}(:,3),{'-'},3,{'r'});
hold on
plot([dates(1) dates(end)],zeros(2,1))
xlim([2008 2013.75])
set(gca,'FontSize',16)
title('Real PCE')

disp('Figure 8 Done')
toc