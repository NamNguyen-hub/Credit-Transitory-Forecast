for ii = 1:3

    lambda = fminsearch(@(lambda0) BN_BVAR_oos_RMSE(y{ii,1},p,lambda0,target_variable(ii)),0,options.optimisation);
        
    [BN_cycle] = BN_BVAR(y{ii,1},p,lambda,target_variable(ii));
    multi_gap(:,ii) = BN_cycle(:,target_variable(ii));
    
end

%% Delete unemployement

y{4,1} = y{2,1}; y{4,1}(:,8) = [];
lambda = fminsearch(@(lambda0) BN_BVAR_oos_RMSE(y{4,1},p,lambda0,target_variable(2)),0,options.optimisation);

[BN_cycle] = BN_BVAR(y{4,1},p,lambda,target_variable(2));
multi_gap(:,4) = BN_cycle(:,target_variable(2));

%%
figure
subplot(2,1,1)
h1=NBERbc(dates,multi_gap(:,1),{':'},3,{'b'});
hold on
h2=plot(dates,multi_gap(:,2),'-.r','LineWidth',3);
hold on
h3=plot(dates,multi_gap(:,3),'-','Color',[.929 .694 .1250],'LineWidth',3);
hold on
plot([dates(1) dates(end)],zeros(2,1));
legend([h1 h2 h3],{'8 Variables','23 Variables','138 Variables'})
set(gca,'FontSize',16)

%%
subplot(2,1,2)
h1=NBERbc(dates,multi_gap(:,2),{'-'},3,{'r'});
hold on
h2=plot(dates,multi_gap(:,4),'-.b','LineWidth',3);
hold on
plot([dates(1) dates(end)],zeros(2,1));
legend([h1 h2],{'Benchmark','Benchmark Omit Umemployment'})
set(gca,'FontSize',16)

disp('Figure 4 Done')
toc
