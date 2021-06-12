function [clean_data] = clearData(modelData, expData)
% Input  : modelData matrix, experimental Data matrix
% Output : Clean modelData if have more data than expData
ncol = size(modelData,2);
N = size(expData,1);
modelTime = modelData(:,1);
clean_data = zeros(N, ncol);

for i = 1:N
    num = expData(i);
    [~, ind] = min(abs(modelTime - num));
    clean_data(i,:) = modelData(ind,:);
end

end