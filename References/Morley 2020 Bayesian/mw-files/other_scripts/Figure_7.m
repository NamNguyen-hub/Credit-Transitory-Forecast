load GDI_data
% just put GDI instead of GDP in the data matrix
GDI_y = [y{2}(:,1) GDI_data  y{2}(:,3:end)];

%%

% run benchmark
lambda = fminsearch(@(lambda0) BN_BVAR_oos_RMSE(y{2},p,lambda0,target_variable(2)),0,options.optimisation);
[BN_cycle] = BN_BVAR(y{2,1},p,lambda,target_variable(2));

% GDI here
lambda = fminsearch(@(lambda0) BN_BVAR_oos_RMSE(GDI_y,p,lambda0,target_variable(2)),0,options.optimisation);
[BN_cycle_tilde] = BN_BVAR(GDI_y,p,lambda,target_variable(2));

%%
figure
h1=NBERbc(dates,BN_cycle(:,2),{'-'},2,{'b'});
hold on
h2=plot(dates,BN_cycle_tilde(:,2),':r','LineWidth',2);
hold on
plot([dates(1) dates(end)],zeros(2,1));
%ylim([-4 4])
set(gca,'FontSize',16)
legend([h1 h2],{'GDP','GDI'},'Location','southeast','Orientation','horizontal')

disp('Figure 7 Done');
toc