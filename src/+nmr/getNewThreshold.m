function [newThr] = getNewThreshold(data,Thr)
%GETNEWTHRESHOLD Summary of this function goes here
%   Detailed explanation goes here
    [X_optimal,pdf,binSize] = nmr.nmr(data);
    N_X = (length(X_optimal)-1)/2;
    newThr = 3;
    for i = 1:N_X
       if(sum(X_optimal(N_X+1+i:end))<1-normcdf(Thr))
           newThr = (i-1)*binSize;
          break; 
       end
    end
%     figure;
%     plot([-N_X:N_X]*binSize,X_optimal,'LineWidth',1.5);
%     hold on;
%     plot([newThr,newThr],[0,max(X_optimal)],'LineWidth',1.5);
end

