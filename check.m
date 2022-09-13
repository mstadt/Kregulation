pars = set_params();
Phi_dtK_ss = 0.08309;

pars.cdKreab_A = 0.294864;%0.000294864*1000; 0.00075*1000; %0.0057; 
pars.cdKreab_B = 0.473015;%0.473015*100; 0.0054*100; %0.0068508;

computed_Phi_cdKreab = 0.04155;
eta_cdKreab = 1; % only changes in MK crosstalk
temp_e = (Phi_dtK_ss-pars.cdKreab_B/100)*pars.cdKreab_A/1000;
temp = 1/(1+exp(temp_e));
cdKreab_ss = (Phi_dtK_ss*temp*eta_cdKreab);

A = 1.1*pars.cdKreab_A;
B = pars.cdKreab_B;
temp2_e = (Phi_dtK_ss-B/100)*A/1000;
temp2 = 1/(1+exp(temp2_e));
cdKreab_NEW = (Phi_dtK_ss*temp2*eta_cdKreab);