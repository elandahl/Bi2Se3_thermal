
%% filmPlotter

clear all;

figure(1);


load thermalFilmOut.m;
L_film=L*1e9; % Convert to nm
 T = [T1' T2'];
Z_tot = [ZZ(:,1)' Z(:,1)'+L]*1e9; % concat,covert to nm
figure(1);clf;hold all; 

ii(1)=1;  % Choose these timepoints
ii(2)=4;
ii(3)=7;
ii(4)=9;
ii(5)=11;
ii(6)=13;

for i = 1:6
  ti = ii(i);
  xlim([0 300]);
  ylim([0 T0*1.15]);

  %legend([num2str(time(ti)*1e9) ' ns'])
  plot(Z_tot,T(ti,:),'LineWidth',3)
end  

% Shade graph, add text
bar([L_film/2],[T0*1.15],L_film/2,'y','LineStyle','none')
xlabel('Depth (nm)','FontSize',16)
ylabel('Temperature (C)','FontSize',16)
text(10,1.07*T0,'Al film','FontSize',16,'FontWeight','bold')
text(10+L_film,1.07*T0,crystal,'FontSize',16,'FontWeight','bold')

% Clever routine for legends
time = time*1e12;
ti=1;
lgd = sprintf('%.0f ps', time(ti));
for idx=2:6, ti = ii(idx); lgd = strvcat(lgd, sprintf('%.0f ps', time(ti))); end,
lgd=cellstr(lgd);
%for idx=1:6, disp(lgd{idx}), end
LEG=legend(lgd);
set(gca, 'FontSize', 16)
set(LEG,'FontSize',16)
%legend([num2str(1e12*time(ii)','%.0f')])
