%% 
clear ; clc ; close all
load data.csv
dataExp =importdata('data.csv');
dataAnane = clearData(dataAn,dataExp);
dataXu2 = clearData(dataXu, dataExp);
dataDewasme = clearData(dataDew, dataExp);

%% Printeo de parametros
% N_param Anane: 15
clc
fprintf("\nData Anane statistics\n", nan);
residual_analysis(dataExp, dataAnane, 15,nan);

% N_param Xu: 13
fprintf("\nData Xu statistics\n",nan)
residual_analysis(dataExp, dataXu2, 13,nan);

% N_param Dewasme: 12
fprintf("\nData dewasme stats\n",nan)
residual_analysis(dataExp, dataDewasme, 12,nan);

%%
time = dataExp(:,1);
data = dataExp(:, 2:end);

subplot(3,1,1)
plot(time, data(:,1), "bo", "MarkerSize", 4)
hold on
plot(time, data(:,2), "ro","MarkerSize", 4)
plot(time, dataAnane(:,2), "LineWidth", 1)
plot(time, dataAnane(:,3), "LineWidth", 1)
hold off

subplot(3,1,2)
plot(time, data(:,3), "bo", "MarkerSize", 4)
hold on
plot(time, dataAnane(:,4), "LineWidth", 1)
hold off
subplot(3,1,3)
plot(time, data(:,4), "bo", "MarkerSize", 4)
hold on
plot(time, dataAnane(:,5), "LineWidth", 1)
hold off