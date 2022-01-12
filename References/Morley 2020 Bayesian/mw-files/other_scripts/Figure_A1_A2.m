% RUN Benchmark again
lambda = fminsearch(@(lambda0) BN_BVAR_oos_RMSE(y{2},p,lambda0,target_variable(2)),0,options.optimisation);
[BN_cycle,Info_decom,shock_decom,FEVD] = BN_BVAR(y{2},p,lambda,target_variable(2),'Decomposition');

%% FEVD (Figure A1)

plot_FEVD = [FEVD.cycle(17,1)./sum(FEVD.cycle(:,1)) FEVD.cycle(1,1)./sum(FEVD.cycle(:,1));...
    FEVD.cycle(17,5)./sum(FEVD.cycle(:,5)) FEVD.cycle(1,5)./sum(FEVD.cycle(:,5))
    FEVD.cycle(17,end)./sum(FEVD.cycle(:,end)) FEVD.cycle(1,end)./sum(FEVD.cycle(:,end))];

subplot(1,2,1)
h1=bar(100*plot_FEVD');
bonecolors = bone(3);
for jj = 1:3
    set(h1(jj),'FaceColor',bonecolors(jj,:),'EdgeColor','k','LineStyle','-');
end
legend({'h = 0','h = 4','h = \infty'},'Location','Northwest')
set(gca,'xticklabel',{'Monetary Policy Shock','Oil Price Shock'})
set(gca,'FontSize',14)
ylim([0 12])
title('Output Gap')
grid on
grid minor

subplot(1,2,2)
bar([FEVD.trend(17) FEVD.trend(1)],0.5,'k');
xlim([0.5 2.5])
set(gca,'xticklabel',{'Monetary Policy Shock','Oil Price Shock'})
set(gca,'FontSize',14)
ylim([0 10])
title('Trend Output')
grid on
grid minor
%% Historical Decomposition 

% Need to mess around with the positive and negative values to stack
figure
NBERbc(dates,BN_cycle(:,target_variable(2)),{'-'},1,{'r'});
J_temp = shock_decom.cycle;
J = J_temp(:,[ 1 17]);
Jneg = J;
Jneg(Jneg>0) = 0;
Jpos = J;
Jpos(Jpos<0) = 0;
h1=bar(dates,Jneg,0.75,'stack');
hold on
h2=bar(dates,Jpos,0.75,'stack');
bonecolors = bone(2);
for jj = 1:2
    set(h1(jj),'FaceColor',bonecolors(jj,:),'EdgeColor','k','LineStyle','-')
    set(h2(jj),'FaceColor',bonecolors(jj,:),'EdgeColor','k','LineStyle','-')
end
hold on
plot(dates,BN_cycle(:,target_variable(2)),'-r','LineWidth',1);
xlim([dates(1) dates(end)])
title('Decomposition of Output Gap into Identified Shocks','FontSize',20)
set(gca,'FontSize',14)
legend(h1,{'Oil Price Shock','Monetary Policy Shock'},'Location','southeast')

disp('Figure A1 & A2 Done')
toc