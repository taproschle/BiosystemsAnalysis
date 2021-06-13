info = readtable('data.csv');
info = table2array(info);

% Se necesita tener Xudata en Workspace
usables = clearData(Xudata, info);
residual_analysis(info(:,2:end),usables(:,2:end), 7);


