shares = std(Info_decom.cycle);

%%
figure
subplot(3,1,1)
bar(1:8,shares(1:8),'w')
hold on
bar(1:8,[0 shares(2) shares(3) 0 0 0 0 shares(8)],'k')
%bar(haves(1:3),shares(haves(1:3)),'r')
set(gca,'xticklabel',Series(1:8),'FontSize',11)
grid on
grid minor
ylim([0 0.5])
xlim([0.5 8.5])

subplot(3,1,2)
bar(1:8,shares(9:16),'w')
hold on
bar(1:8,[0 shares(10) 0 0 shares(13) 0 0 0],'k')
set(gca,'xticklabel',Series(9:16),'FontSize',11)
grid on
grid minor
ylim([0 0.5])
xlim([0.5 8.5])

subplot(3,1,3)
bar(1:7,shares(17:23),'w')
hold on
bar(1:7,[shares(17) 0 shares(19) 0 0 0 shares(23)],'k')
set(gca,'xticklabel',Series(17:23),'FontSize',11)
grid on
grid minor
ylim([0 0.5])
xlim([0.5 7.5])


disp('Figure 3 Done')
toc
