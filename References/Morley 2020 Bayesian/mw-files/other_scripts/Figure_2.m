lambda = fminsearch(@(lambda0) BN_BVAR_oos_RMSE(y{2},p,lambda0,target_variable(2)),0,options.optimisation);
[BN_cycle,Info_decom,shock_decom,FEVD,Std_err] = BN_BVAR(y{2},p,lambda,target_variable(2),'Decomposition','Standard_Errors');


%%
figure
h1=NBERbc(dates,BN_cycle(:,target_variable(2)),{'-'},3,{'r'});
hold on
plot(dates,BN_cycle(:,target_variable(2))+1.645*Std_err(target_variable(2)),'--k','LineWidth',1);
hold on
plot(dates,BN_cycle(:,target_variable(2))-1.645*Std_err(target_variable(2)),'--k','LineWidth',1);
hold on
plot([dates(1) dates(end)],zeros(2,1),'-k','LineWidth',2)
ylim([-8 6])
set(gca,'FontSize',16)

disp('Figure 2 Done')
toc