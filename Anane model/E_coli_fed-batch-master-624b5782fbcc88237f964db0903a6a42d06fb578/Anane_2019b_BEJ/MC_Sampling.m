function MonteCarlo_data    = MC_Sampling(N,Measurements,Standard_devs) % correlated parameter samples

%E Anane
% TU-Berlin, 20th June 2018

%% Check for errors in data
if length(Measurements)~=length(Standard_devs)
     error('Each measurement must have a standard deviation')
 end
  
 zero_indices = find(~Standard_devs);
 if any(Measurements(zero_indices))~=0; %#ok<FNDSB>
     error('A non-zero measurement should not have a 0 standard deviation')   
 end 
 
%% make normal distribution from measurements and their std deviations
pd                 = cell(1,length(Measurements));

    for j1         = 1:length(Measurements)
        pd{:,j1}   = makedist('Normal','mu',Measurements(j1),'sigma',Standard_devs(j1));  %Generate probability dist. functions with measurements as means and standard deviation from triplicate samples. 
    end
    
%% Generate random samples from the normal distribution spaces
l = length(pd);

x = lhsdesign(N,l);                         % generate latin hypercube samples. See MATLAB documentation for more information)
Generated_data=zeros(N,l);                  % preallocation for the matrix

    for i=1:l           
        prob    =   x(:,i);
        Generated_data(:,i) = icdf('Normal',prob,Measurements(i),Standard_devs(i));   % map latin hypercube samples to real values using inverse cumulative distribution functions
    end

%% Replace NaNs (for data points with standard deviation = 0) with 0s.
Generated_data(isnan(Generated_data)) = 0;

MonteCarlo_data = Generated_data;

end


