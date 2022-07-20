function F0 = getDFmovmin(datIn,window)
    F0 = datIn;
    repeatTimes = 10;
    for i = 1:repeatTimes
        F0 = min(datIn,F0);
        F0 = movmean(F0,window,3);
    end
    
%     dF = datIn-F0;
end

