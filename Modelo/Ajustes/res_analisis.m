function [stistics] = res_analisis(model)


% To know the 
kL = k.*0.9999999999;
kU = k.*1.0000000001 + ones(1,length(k))*1e-15;
[k2,~,res,~,~,~,Jac] = lsqcurvefit(@lsq_func,k,texp,yexp,kL,kU);
diffs = sum((k-k2).^2);

% AIC & AICc calculations
SS    = sum(sum(res.^2)); %First col and then rows
m     = length(k2);
[N,~] = size(res);
AIC   = N*log(SS/N)+2*(m+1);
AICc  = AIC + 2*(m+1)*(m+2)/(N-m-2);

end