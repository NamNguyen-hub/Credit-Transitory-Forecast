lambda_grid = [1e-20 0.01:0.01:0.19 0.4:0.01:0.55 0.6:0.1:1]';

for kk = 1:3    
    for jj = 1:size(lambda_grid,1)
        oos_RMSE(jj,kk) = BN_BVAR_oos_RMSE(y{kk,1},p,lambda_grid(jj),target_variable(kk));
    end    
end

AR1_RMSE = oos_forecast_error_ar(y{1}(:,target_variable(1)),1);

%%
figure
plot(lambda_grid,oos_RMSE(:,1),':','LineWidth',3);
hold on
plot(lambda_grid,oos_RMSE(:,2),'--','LineWidth',3);
hold on
plot(lambda_grid,oos_RMSE(:,3),'-','LineWidth',3);
hold on
plot([lambda_grid(1) lambda_grid(end)],[AR1_RMSE AR1_RMSE],'-k');
xlabel('\lambda','FontSize',16)
legend('8 Variables','23 Variables','138 Variables')
set(gca,'FontSize',16)

disp('Figure 5 Done')
toc
