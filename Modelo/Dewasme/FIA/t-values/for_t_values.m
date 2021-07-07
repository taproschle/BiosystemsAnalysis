% t value multiple times calculated
clear ; clc ; close all

times = 100;
adj_params = 2;
ts = ones(times,adj_params);

for i = 1:times
    run Main.m
    for j = 1:adj_params
        ts(i,j) = t_values(j);
    end
end

table = ones(adj_params,2);
for i = 1:adj_params
    eval(sprintf('k%dm = mean(ts(:,%d));',i,i))
    eval(sprintf('table(%d,1) = k%dm;',i,i))
    eval(sprintf('k%dstd = std(ts(:,%d));',i,i))
    eval(sprintf('table(%d,2) = k%dstd;',i,i))    
end
