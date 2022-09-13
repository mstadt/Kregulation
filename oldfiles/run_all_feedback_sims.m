% doFF, doINS, doALDNKA, doALDsec
sims = [1,1,1,1; % all feedbacks on 1
        0,1,1,1; % GI FF off 2
        1,0,1,1; % insulin off 3
        1,1,0,1; % ALD NKA off 4
        1,1,1,0; % ALD sec off 5
        0,0,0,0; % all off 6
        1,0,0,0; % only GI FF 7
        0,1,0,0; % only insulin 8
        0,0,1,0; % only ALD NKA 9 
        0,0,0,1; % only ALD sec 10
        ];

days = 50;

fname = strcat('./results/', date, '_feedback_simulations');
X = {};
T = {};
for ii = 1:size(sims,1)
    fprintf('simulation %i \n', ii)
    [T{ii}, X{ii}] = feedback_sim(sims(ii, :), days);
end
save(fname, 'T', 'X', 'days', 'sims')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T, X] = feedback_sim(sim, days)
    doFF     = sim(1);
    doINS    = sim(2);
    doALDNKA = sim(3);
    doALDsec = sim(4);


    pars = set_params();
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