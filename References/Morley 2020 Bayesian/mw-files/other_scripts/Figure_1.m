
fig1_dates = (1959.25:0.25:2016.75)';
rawdata = xlsread('small_var_data.xlsx','B7:I238');

x = 100*diff(log(rawdata));
x = [x(:,1) rawdata(2:end,2:3) x(:,4:8)];


%%

BN_cycle_temp1 = multivariate_BN_ols(x(:,1),2);
BN_cycle_temp2 = multivariate_BN_ols(x(:,1:2),p);
BN_cycle_temp3 = multivariate_BN_ols(x(:,[1:2 4:5]),p);
BN_cycle_temp4 = multivariate_BN_ols(x,p);

%%
figure
h1=NBERbc(fig1_dates,BN_cycle_temp1(:,1),{'-'},3,{'r'});
hold on
h2=plot(fig1_dates,BN_cycle_temp2(:,1),'-.','LineWidth',3);
hold on
h3=plot(fig1_dates,BN_cycle_temp3(:,1),':','LineWidth',3);
hold on
plot([fig1_dates(1) fig1_dates(end)],zeros(2,1),'-k','LineWidth',2)
ylim([-6 5])
legend([h1 h2 h3],{'Univariate','Two Variable VAR','Four Variable VAR'})
set(gca,'FontSize',16)

disp('Figure 1 Done')
toc

