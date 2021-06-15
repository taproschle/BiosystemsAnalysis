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
outAnane = residual_analysis(info(:,2:end),Aninfo(:,2:end), 15, 'Anane');
outXu = residual_analysis(info(:,2:end),Xuinfo(:,2:end), 13, 'Xu');
outDew = residual_analysis(info(:,2:end),Dewinfo(:,2:end), 12, 'Dewasme');

%% 
resAnane = outAnane.res;
resXu = outXu.res;
resDew = outDew.res;

%% 
clc
cols = size(resAnane,2);
pvalues = zeros(cols, 3);
isNorm = zeros(cols, 3);
orden = ['X', 'E', 'S', 'O'];
for i = 1:4
    [rAn,pAn] = adtest(resAnane(:,i));
    [rXu,pXu] = adtest(resXu(:,i));
    [rDew,pDew] = adtest(resDew(:,i));
    pvalues(i,:) = [pAn pXu pDew];
    isNorm(i,:) = [rAn rXu rDew];
end

%%

plot(resAnane(:,1), "bo");
%plot(resXu(:,1), "ro");
hold on
yline(0, 'k-', "LineWidth", 0.7)
hold off 






