function riseT = getRisingTime(x0,t0,t1,thrVec)
    [maxV,tMax] = max(x0(t0:t1));
    tMax = tMax + t0 - 1;   
    
    [~,tMin] = min(x0(t0:tMax));
    t0 = t0 + tMin - 1;
    % extenstion
    t0 = max(t0 - 10,1);
    [minV,tMin] = min(x0(t0:tMax));
    tMin = tMin + t0 - 1;
    
    if(tMin==tMax)
        riseT = tMax;
        return;
    end

%     figure;plot(x0);hold on;plot(tMin:t1,x0(tMin:t1));
    
    x0 = x0(tMin:tMax);
    x0 = (x0-minV)/(maxV-minV);

    try
        riseTs = zeros(numel(thrVec),1);
        for i = 1:numel(thrVec)
            tAch = find(x0>thrVec(i),1);
            riseTs(i) = tAch;
        end
        riseT = mean(riseTs) + tMin-1;
    catch
        riseT = tMin;
    end
end