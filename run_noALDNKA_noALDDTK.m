% run simulation experiements with and without FF effect
clear all;
close all;
% IDEA: run for 1 day of 3 meals....
%           or maybe a week of normal intake?

days = 14; % number of days to run the simulation for
%%%%%%%%%%%%%%%%%%
% simulation 1
% baseline model
%%%%%%%%%%%%%%%%%%
disp('**sim 1**')
pars1 = set_params();

Kin1.Kin_type = 'long_simulation'; % TO DO: CHECK WHAT THIS DOES!
Kin1.Meal = 1;
Kin1.KCL = 1;

MealInfo1.t_breakfast = 7;
MealInfo1.t_lunch = 13;
MealInfo1.t_dinner = 19;
MealInfo1.K_amount = 100/3; % how much K is ingested PER MEAL (3 times per day)
MealInfo1.meal_type = 'normal'; % TO DO: change to a "normal" meal type

do_FF1 = 1;
do_ins1 = 1;

% sim 1 SS condition
disp('get sim 1 SS')
IG_file1 = './IGdata/KregSS.mat';
[SSdata1, exitflag1, residual1] = getSS(IG_file1, pars1, Kin1, ...
                                           'MealInfo', {MealInfo1.t_breakfast, MealInfo1.t_lunch, MealInfo1.t_dinner,...
                                                        MealInfo1.K_amount, MealInfo1.meal_type}, ...
                                            'do_insulin', [do_ins1, pars1.insulin_A, pars1.insulin_B], ...
                                            'do_FF', [do_FF1, pars1.FF]);
if exitflag1 <=0
    disp('residuals')
    disp(residual1)
    error('******sim 1 SS did not converge*****')
end
disp('sim 1 SS finished')

%  run simulation 1
opts = odeset('MaxStep', 20);
x0 = SSdata1;
x_p0 = zeros(size(SSdata1));
t0 = 0;
tf = days*1440; 
tspan = t0:0.5:tf;

disp('start sim 1')
[T1,X1] = ode15i(@(t,x,x_p) k_reg_mod(t,x,x_p, pars1, ...
                                'Kin_type', {Kin1.Kin_type, Kin1.Meal, Kin1.KCL}, ...
                                'MealInfo', {MealInfo1.t_breakfast, MealInfo1.t_lunch, MealInfo1.t_dinner,...
                                                        MealInfo1.K_amount, MealInfo1.meal_type}, ...
                                'do_insulin', [do_ins1, pars1.insulin_A, pars1.insulin_B], ...
                                'do_FF', [do_FF1, pars1.FF]), ...
                            tspan, x0, x_p0, opts);
disp('sim 1 finished')

%%%%%%%%%%%%%%%%%%
% simulation 2
% ALD effect off
%%%%%%%%%%%%%%%%%%
disp('**sim 2**')
pars2 = set_params();

Kin2.Kin_type = 'long_simulation'; % TO DO: CHECK WHAT THIS DOES!
Kin2.Meal = 1;
Kin2.KCL = 1;

MealInfo2.t_breakfast = 7;
MealInfo2.t_lunch = 13;
MealInfo2.t_dinner = 19;
MealInfo2.K_amount = 100/3; % how much K is ingested PER MEAL (3 times per day)
MealInfo2.meal_type = 'normal'; % TO DO: change to a "normal" meal type

do_ins2 = 1;
do_FF2 = 1; % 
do_ALD_NKA2 = 0; % ALD effect on NKATPase
do_ALD_sec2 = 1; % ALD effect on DT K sec

% sim 2 SS condition

disp('get sim 2 SS')
IG_file2 = './IGdata/KregSS.mat';
[SSdata2, exitflag2, residual2] = getSS(IG_file2, pars2, Kin2, ...
                                           'MealInfo', {MealInfo2.t_breakfast, MealInfo2.t_lunch, MealInfo2.t_dinner,...
                                                        MealInfo2.K_amount, MealInfo2.meal_type},...
                                            'do_insulin', [do_ins2, pars2.insulin_A, pars2.insulin_B], ...
                                            'do_FF', [do_FF2, pars2.FF], ...
                                            'do_ALD_NKA', do_ALD_NKA2, ...
                                            'do_ALD_sec', do_ALD_sec2);
if exitflag2 <=0
    disp('residuals')
    disp(residual2)
    error('******sim 2 SS did not converge*******')
end
disp('sim 2 SS finished')

%  run simulation 2
opts = odeset('MaxStep', 20);
x0 = SSdata2;
x_p0 = zeros(size(SSdata2));
t0 = 0;
tf = days*1440; 
tspan = t0:0.5:tf;

disp('start sim 2')
[T2,X2] = ode15i(@(t,x,x_p) k_reg_mod(t,x,x_p, pars2, ...
                                'Kin_type', {Kin2.Kin_type, Kin2.Meal, Kin2.KCL}, ...
                                'MealInfo', {MealInfo2.t_breakfast, MealInfo2.t_lunch, MealInfo2.t_dinner,...
                                                        MealInfo2.K_amount, MealInfo2.meal_type}, ...
                                'do_insulin', [do_ins2, pars2.insulin_A, pars2.insulin_B], ...
                                'do_FF', [do_FF2, pars2.FF], ...
                                'do_ALD_NKA', do_ALD_NKA2, ...
                                'do_ALD_sec', do_ALD_sec2), ...
                            tspan, x0, x_p0, opts);
disp('sim 2 finished')

%%%%%%%%%%%%%%%%%%
% simulation 3
% both ALD effects off
%%%%%%%%%%%%%%%%%%
disp('**sim 3**')
pars3 = set_params();

Kin3.Kin_type = 'long_simulation'; % TO DO: CHECK WHAT THIS DOES!
Kin3.Meal = 1;
Kin3.KCL = 1;

MealInfo3.t_breakfast = 7;
MealInfo3.t_lunch = 13;
MealInfo3.t_dinner = 19;
MealInfo3.K_amount = 100/3; % how much K is ingested PER MEAL (3 times per day)
MealInfo3.meal_type = 'normal'; % TO DO: change to a "normal" meal type

do_ins3 = 1; % no insulin effect
do_FF3 = 1; % no FF effect
do_ALD_NKA3 = 0; % ALD effect on NKATPase
do_ALD_sec3 = 0; % ALD effect on DT K sec

% sim 3 SS condition

disp('get sim 3 SS')
IG_file3 = './IGdata/KregSS.mat';
[SSdata3, exitflag3, residual3] = getSS(IG_file3, pars3, Kin3, ...
                                           'MealInfo', {MealInfo3.t_breakfast, MealInfo3.t_lunch, MealInfo3.t_dinner,...
                                                        MealInfo3.K_amount, MealInfo3.meal_type},...
                                            'do_insulin', [do_ins3, pars3.insulin_A, pars3.insulin_B], ...
                                            'do_FF', [do_FF3, pars3.FF], ...
                                            'do_ALD_NKA', do_ALD_NKA3, ...
                                            'do_ALD_sec', do_ALD_sec3);
if exitflag3 <=0
    disp('residuals')
    disp(residual3)
    error('******sim 3 SS did not converge*******')
end
disp('sim 3 SS finished')

%  run simulation 3
opts = odeset('MaxStep', 20);
x0 = SSdata3;
x_p0 = zeros(size(SSdata3));
t0 = 0;
tf = days*1440; 
tspan = t0:0.5:tf;

disp('start sim 3')
[T3,X3] = ode15i(@(t,x,x_p) k_reg_mod(t,x,x_p, pars3, ...
                                'Kin_type', {Kin3.Kin_type, Kin3.Meal, Kin3.KCL}, ...
                                'MealInfo', {MealInfo3.t_breakfast, MealInfo3.t_lunch, MealInfo3.t_dinner,...
                                                        MealInfo3.K_amount, MealInfo3.meal_type}, ...
                                'do_insulin', [do_ins3, pars3.insulin_A, pars3.insulin_B], ...
                                'do_FF', [do_FF3, pars3.FF], ...
                                'do_ALD_NKA', do_ALD_NKA3, ...
                                'do_ALD_sec', do_ALD_sec3), ...
                            tspan, x0, x_p0, opts);
disp('sim 3 finished')


% plot simulation
do_plt = 1;
if do_plt
    clear T X params Kin_opts MealInfo
    disp('**plotting simulation results**')
    Kin_opts{1} = Kin1;
    Kin_opts{2} = Kin2;
    Kin_opts{3} = Kin3;
    MealInfo{1} = MealInfo1;
    MealInfo{2} = MealInfo2;
    MealInfo{3} = MealInfo3;
    params{1} = pars1;
    params{2} = pars2;
    params{3} = pars3;
    T{1} = T1;
    T{2} = T2;
    T{3} = T3;
    X{1} = X1;
    X{2} = X2;
    X{3} = X3;
    labels{1} = 'baseline model';
    labels{2} = 'no ALD NKA effect';
    labels{3} = 'no ALD effects';
    plot_noALD2_sim(T,X,params,Kin_opts,labels,tf,MealInfo,days)
end

