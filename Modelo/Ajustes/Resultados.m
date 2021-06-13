%% 

dataAnane = clearData(dataAn,dataExp);
dataXu2 = clearData(dataXu, dataExp);
dataDewasme = clearData(dataDew, dataExp);

%% Printeo de parametros
% N_param Anane: 15
clc
fprintf("\nData Anane statistics\n")
residual_analysis(dataExp, dataAnane, 15);

% N_param Xu: 13
fprintf("\nData Xu statistics\n")
residual_analysis(dataExp, dataXu2, 13);

% N_param Dewasme: 12
fprintf("\nData dewasme stats\n")
residual_analysis(dataExp, dataDewasme, 12);