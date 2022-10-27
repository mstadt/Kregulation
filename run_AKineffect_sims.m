% run simulation experiements with and without FF effect
clear all;
%%%%    Begin User Input %%%%%%%
days = 10; % number of days to run the simulation for
% sim 1 options
label1 = 'baseline model';
AKin1_per  = 1.0;
% sim 2 options
label2 = 'A_{Kin} 5% decrease';
AKin2_per  = 0.95;
% sim 3 options
label3 = 'A_{Kin} 10% decrease';
AKin3_per = 0.9;
% sim 4 options
label4 = 'A_{Kin} 5% increase';
AKin4_per  = 1.05;
% sim 5 options
label5 = 'A_{Kin} 10% increase';
AKin5_per = 1.1;
%%%% end user input %%%%%%%%

AKin_sims = [AKin1_per;
        AKin2_per;
        AKin3_per;
        AKin4_per;
        AKin5_per];


X = {};
T = {};
days = 10;

for ii = 1:size(AKin_sims,1)
    fprintf('simulation %i \n', ii)
    [T{ii}, X{ii}] = AKin_sim(AKin_sims(ii,:), days);
end

labels = {label1, label2,label3, label4,label5};
MealInfo.t_breakfast = 7;
MealInfo.t_lunch = 13;
MealInfo.t_dinner = 19;
MealInfo.K_amount = 100/3;
MealInfo.meal_type = 'normal';
Kin.Kin_type = 'long_simulation';
Kin.Meal = 1;
Kin.KCL = 1;
pars = set_params();

plot_5_sims(T,X,pars,Kin,labels,MealInfo,days)


function [T, X] = AKin_sim(AKin_perchange, days)
    pars = set_params();
    pars.FF = AKin_perchange*pars.FF; % change A_Kin parameter

    Kin.Kin_type = 'long_simulation';
    Kin.Meal = 1;
    Kin.KCL = 1;

    MealInfo.t_breakfast = 7;
    MealInfo.t_lunch = 13;
    MealInfo.t_dinner = 19;
    MealInfo.K_amount = 100/3;
    MealInfo.meal_type = 'normal';
    disp('get SS')
    IG_file = './IGdata/KregSS.mat';

    doINS = 1;
    doFF = 1;
    doALDNKA = 1;
    doALDsec = 1;
    [SSdata, exitflag, residual] = getSS(IG_file, pars, Kin, ...
                                            'MealInfo', {MealInfo.t_breakfast, MealInfo.t_lunch, MealInfo.t_dinner, ...
                                                            MealInfo.K_amount, MealInfo.meal_type}, ...
                                            'do_insulin', [doINS, pars.insulin_A, pars.insulin_B], ...
                                            'do_FF', [doFF, pars.FF], ...
                                            'do_ALD_NKA', doALDNKA, ...
                                            'do_ALD_sec', doALDsec);
    if exitflag <= 0
        disp('residuals')
        disp(residual)
        error('sim SS did not converge')
    end
    disp('SS solved')

    % run simulation
    opts = odeset('MaxStep', 20);
    x0 = SSdata;
    x_p0 = zeros(size(SSdata));
    t0 = 0;
    tf = days*1440;
    tspan = t0:0.5:tf;
    disp('run sim')
    [T, X] = ode15i(@(t,x,x_p) k_reg_mod(t,x,x_p,pars,...
                                            'Kin_type', {Kin.Kin_type, Kin.Meal, Kin.KCL}, ...
                                            'MealInfo', {MealInfo.t_breakfast, MealInfo.t_lunch, MealInfo.t_dinner, ...
                                                                MealInfo.K_amount, MealInfo.meal_type}, ...
                                            'do_insulin', [doINS, pars.insulin_A, pars.insulin_B], ...
                                            'do_FF', [doFF, pars.FF], ...
                                            'do_ALD_NKA', doALDNKA, ...
                                            'do_ALD_sec', doALDsec), ...
                                        tspan, x0, x_p0, opts);
    disp('simulation finished')
end % feedback sim



