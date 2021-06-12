function out = residual_analysis(expdata, modeldata, num_param, show)

if nargin < 4; show = true; end
if nargin < 3; error('Insuficientes parametros'); end

[NExp, nExp] = size(modeldata);
[Ndata, ndata] = size(expdata);

if NExp ~= Ndata || nExp < ndata; error('Deben ser del mismo'); end

modeldata = modeldata(:,1:ndata);

RSS = (sum(sum(abs(expdata - modeldata))))^2;
means = mean(expdata);

TSS = sum(sum(abs(expdata - repmat(means, Ndata, 1))))^2;

out.RSq    = 1 - RSS/TSS;
out.RSqAdj = 1 - RSS/TSS*(NExp - num_param - 1)/(NExp - 1);

out.AIC   = NExp*log(RSS/NExp)+2*(num_param+1);
out.AICc  = out.AIC + 2*(num_param+1)*(num_param+2)/(NExp-num_param-2);

out.BIC = NExp*log(RSS/NExp) + num_param*log(NExp);

if show == true
   fprintf('Statisctic of the model\n')
   fprintf('Nparameters : %f \n', num_param )
   fprintf('Rsquared    : %f \n', out.RSq)
   fprintf('Rsquared Adj: %f \n', out.RSqAdj)
   fprintf('AIC         : %f \n', out.AIC)
   fprintf('AICc        : %f \n', out.AICc)
   fprintf('BIC         : %f \n', out.BIC)
end

end