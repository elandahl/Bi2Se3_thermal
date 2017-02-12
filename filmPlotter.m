
%% filmPlotter

clear all;
L = 20e-9; % film thickness in m
fluence = 0.1; % fluence in mJ/cm^2
times = (1e-10:2e-10:1e-8); % in seconds
[T1 T2 z1 z2] = Bi2Se3_thermal (L, fluence, times, 1e-5);

figure(1);
clear all;
load thermalFilmOut.m; % If not already called
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
%text(10,1.07*T0,'Film','FontSize',16,'FontWeight','bold')
%text(10+L_film,1.07*T0,'Substrate','FontSize',16,'FontWeight','bold')

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
hold off;

figure(2);clf;hold on
for i = 1:length(time)
  T1_avg(i) = mean(T1(:,i));
end
alpha_t = 1.9e-5; % Thermal expansion, 1/K
dTheta = -(180/pi)*T1_avg*alpha_t/cot(7.7*pi/180); %centroid shift in deg
plot(time*1e-3,1000*dTheta,'LineWidth',3)
xlabel('Time (ns)','FontSize',16)
ylabel('Peak centroid shift (mdeg)','FontSize',16)
set(gca, 'FontSize', 16)
hold off;
 
