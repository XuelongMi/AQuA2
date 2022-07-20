function [F0] = baselineLinearEstimate(datIn,cut,movAvgWin,nonvalid,ff)
    
    % non-valid is to show the region caused by registration edges
    % grow 5 time points for nan valid part in temporal
    [H,W,T] = size(datIn);
    datMA = datIn;
    datMA(nonvalid) = nan;
    datMA = movmean(datMA,movAvgWin,3,'omitnan');
    datMA(nonvalid) = nan;
    datMA = reshape(datMA,[],T);
    step = round(0.5*cut);
    maxV = max(datIn(:));

    nSegment = ceil(T/step)-1;
    minPosition = zeros(H*W,nSegment);
    for k = 1:nSegment
        if(exist('ff','var'))
            waitbar(0.25*k/nSegment,ff,'Detecting active regions...');
        end
        t0 = 1 + (k-1)*step;
        t1 = min(T,t0+cut);
        [~,curP] = nanmin(datMA(:,t0:t1),[],2);
        minPosition(:,k) = curP+t0-1;
    end
    
    F0 = zeros(size(datMA),'single');
    parfor i = 1:H*W
        curP = unique(minPosition(i,:));
        value = datMA(i,curP);
        curP = curP(~isnan(value));
        value = value(~isnan(value));
        nMin = numel(value);
        curve = zeros(1,T);
        if(nMin==0)
            % all nonValid
            curve = maxV;
        else
           % first part
            curve(1:curP(1)) = value(1);
            % end part
            curve(curP(nMin):T) = value(nMin);
            % middle part
            for k = 1:nMin-1
                mt1 = curP(k);
                mt2 = curP(k+1);
                curve(mt1:mt2) = value(k) + (value(k+1)-value(k))/(mt2-mt1)*[0:mt2-mt1]; 
            end
        end
        F0(i,:) = curve;
    end
    
    F0(nonvalid) = maxV;
    F0 = reshape(F0,[H,W,T]);
end