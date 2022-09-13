function plot_ins_sim(T,X,params,Kin_opts,labels,tf,MealInfo,days)
%close all
t1 = T{1};
t2 = T{2};

    % color options
c1= [0.4940 0.1840 0.5560]; % purple [0.3010 0.7450 0.9330]; % blue 
c2 = [0.4660 0.6740 0.1880]; %green[0.4940 0.1840 0.5560] ; % green %[0.6350 0.0780 0.1840]; % maroon[0.3010 0.7450 0.9330]; 
c3 = [0.3010 0.7450 0.9330]; % blue
c4 = [0.4660 0.6740 0.1880]; %green %[0.4940 0.1840 0.5560]; % purple
lightGrey2   = [0.7 0.7 0.7];

% fontsizes
f.title = 16;
f.xlabel = 12;
f.ylabel = 12;
f.legend = 14;
f.ticks = 16;
f.lab = 12;
mark_size1 = 7;
mark_size2 = 20;
lw = 1.5; % linewidth
% markers
m1 = 'o';
m2 = '^';
m3 = 'd';
m4 = 's';
% linestyles
ls1 = '-';
ls2 = '-';

labx = 0.25;
laby = 0.95;
% figure
xmin = 0;
xmax = 14;
figure
% Phi_Kin
s = subplot(3,2,1);
Kinvals1 = zeros(size(T{1}));
for ii = 1:length(T{1})
    [Kinvals1(ii), ~] = get_PhiKin(T{1}(ii), 0, params{1}, Kin_opts{1}, MealInfo{1});
end % ii
plot(s, t1/1440, Kinvals1, 'linewidth', lw, 'color', 'k')
hold on
xlabel('time (days)', 'fontsize', f.xlabel)
ylabel('mEq/min', 'fontsize', f.ylabel)
title('K^+ intake (\Phi_{Kin})', 'fontsize', f.title)
xlim([xmin, xmax])
ylim([-0.1,2])
text(xmin + labx,gca().YLim(1) + laby*(gca().YLim(2) - gca().YLim(1)), '(a)', 'fontsize', f.lab)
hold off

% K_plasma
s = subplot(3,2,3);
varnum = 5;
vals1 = X{1}(:,varnum);
vals2 = X{2}(:,varnum);
plot(s,t1/1440, vals1, 'linewidth', lw, 'color', c1, 'linestyle', ls1)
hold on
plot(s,t2/1440, vals2, 'linewidth', lw, 'color', c2, 'linestyle', ls2)
yline(3.5, 'linewidth', 1.5, 'color', lightGrey2, 'linestyle', '-')
yline(5.0, 'linewidth', 1.5, 'color', lightGrey2, 'linestyle', '-')
xlabel('time (days)', 'fontsize', f.xlabel)
ylabel('mEq/L', 'fontsize', f.ylabel)
title('Plasma [K^+] (K_{plasma})', 'fontsize', f.title)
xlim([xmin, xmax])
ylim([3.5,5.3])
text(xmin + labx,gca().YLim(1) + laby*(gca().YLim(2) - gca().YLim(1)), '(c)', 'fontsize', f.lab)
hold off

% K_IC
s = subplot(3,2,5);
varnum = 8;
vals1 = X{1}(:,varnum);
vals2 = X{2}(:,varnum);
plot(s,t1/1440, vals1, 'linewidth', lw, 'color', c1, 'linestyle', ls1)
hold on
plot(s,t2/1440, vals2, 'linewidth', lw, 'color', c2, 'linestyle', ls2)
yline(120, 'linewidth', 1.5, 'color', lightGrey2, 'linestyle', '-')
yline(140, 'linewidth', 1.5, 'color', lightGrey2, 'linestyle', '-')
xlabel('time (days)', 'fontsize', f.xlabel)
ylabel('mEq/L', 'fontsize', f.ylabel)
title('Intracellular [K^+] (K_{IC})', 'fontsize', f.title)
xlim([xmin, xmax])
ylim([120, 150])
text(xmin + labx,gca().YLim(1) + laby*(gca().YLim(2) - gca().YLim(1)), '(e)', 'fontsize', f.lab)
hold off

% Phi_dtsec
s = subplot(3,2,2);
varnum = 18;
vals1 = X{1}(:,varnum);
vals2 = X{2}(:,varnum);
plot(s, t1/1440, vals1, 'linewidth', lw, 'color', c1, 'linestyle', ls1)
hold on
plot(s, t2/1440, vals2, 'linewidth', lw, 'color', c2, 'linestyle', ls2)
xlabel('time (days)', 'fontsize', f.xlabel)
ylabel('mEq/min', 'fontsize', f.ylabel)
title('DT K^+ secretion (\Phi_{dt-Ksec})', 'fontsize', f.title)
xlim([xmin, xmax])
ylim([0, 0.32])
text(xmin + labx,gca().YLim(1) + laby*(gca().YLim(2) - gca().YLim(1)), '(b)', 'fontsize', f.lab)
hold off

% Phi_cd transport
s = subplot(3,2,4);
varnum1 = 23;
varnum2 = 26;
vals1 = X{1}(:,varnum1) - X{1}(:,varnum2);
vals2 = X{2}(:,varnum1) - X{2}(:,varnum2);
plot(s, t1/1440, vals1, 'linewidth', lw, 'color', c1, 'linestyle', ls1)
hold on
plot(s, t2/1440, vals2, 'linewidth', lw, 'color', c2, 'linestyle', ls2)
xlabel('time (days)', 'fontsize', f.xlabel)
ylabel('mEq/min', 'fontsize', f.ylabel)
title('CD K^+ transport (\Phi_{cd-Ksec} - \Phi_{cd-Kreab})', 'fontsize', f.title)
xlim([xmin, xmax])
ylim([-0.20,0.0])
text(xmin + labx,gca().YLim(1) + laby*(gca().YLim(2) - gca().YLim(1)), '(d)', 'fontsize', f.lab)
hold off

% Phi_uK
s = subplot(3,2,6);
varnum = 28;
vals1 = X{1}(:,varnum);
vals2 = X{2}(:,varnum);
plot(s, t1/1440, vals1, 'linewidth', lw, 'color', c1, 'linestyle', ls1)
hold on
plot(s, t2/1440, vals2, 'linewidth', lw, 'color', c2, 'linestyle', ls2)
xlabel('time (days)', 'fontsize', f.xlabel)
ylabel('mEq/min', 'fontsize', f.ylabel)
title('Urinary K^+ excretion (\Phi_{uK})', 'fontsize', f.title)
xlim([xmin, xmax])
ylim([0,0.2])
text(xmin + labx,gca().YLim(1) + laby*(gca().YLim(2) - gca().YLim(1)), '(f)', 'fontsize', f.lab)
hold off
    
    

legend(labels{1}, labels{2}, 'fontsize', f.legend)
end % plot_FF_sim