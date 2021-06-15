% read of the data
info = readtable('data.csv');
info = table2array(info);

load simdata_xu
load simdata_dew
load simdata_an

%%
Xuinfo = clearData(dataXu, info);
Dewinfo = clearData(dataDew, info);
Aninfo = clearData(dataAn, info);

%% Analisis
clc
residual_analysis(info(:,2:end),Aninfo(:,2:end), 15, 'Anane');
residual_analysis(info(:,2:end),Xuinfo(:,2:end), 13, 'Xu');
residual_analysis(info(:,2:end),Dewinfo(:,2:end), 12, 'Dewasme');
