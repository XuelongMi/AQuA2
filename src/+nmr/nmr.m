function [X_optimal,pdf,binSize] = nmr(data,weightType,nonincreasing)
% Recover noise model from the data
if(~exist('weightType','var'))
    weightType = 'Average';
end
if(~exist('nonincreasing','var'))
    nonincreasing = 0;
end

% The last dimension denotes time dimension
sz = size(data);
if(length(sz)==2&& sz(2)==1)    % for column vector
    dif = data(2:end)-data(1:end-1);
else
    data = reshape(data,[],sz(end));
    dif = data(:,2:end) - data(:,1:end-1);
end
[pdf,binSize] = nmr.getDistribution(dif);
X_optimal = nmr.pdf_recover(pdf,numel(data(:)),weightType,nonincreasing);

end
