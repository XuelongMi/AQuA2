function [mergingInfo,majorityEvt] = createMergingInfo(evtLst,majorityEvt,ccRegions,dF,opts)
    [H,W,T] = size(dF);
    dF = reshape(dF,[],T);
    N = numel(evtLst);
    overlap = opts.overlap;
    gap = 2*opts.smoT+1;
    %% calculate peak times
    for i = 1:N
        ihw = majorityEvt{i}.ihw;
        TW = majorityEvt{i}.TW;
        curve = mean(dF(ihw,:),1);
        s0 = sqrt(median((curve(2:end)-curve(1:end-1)).^2)/0.9133);
        curve = curve(TW);
        % new method
%         [maxV,tPeak] = max(curve);
%         [minV,tMin] = min(curve(1:tPeak));
%         % 50% rising
%         riseT = find(curve>(maxV-minV)*0.5 + minV,1) + TW(1)-1;
%         if(isempty(riseT))
%             riseT = TW(1);
%         end
        
        smo = 1;
        if(numel(TW)>1)
            while(1)
                if(smo>0)
                    curveSmo = imgaussfilt(curve,smo);
                else
                    curveSmo = curve;
                end
                [maxV,~] = max(curveSmo);
    %                 tPeak = tPeak + t0 - 1;
                left = false(numel(TW),1);
                left(1) = true;
                left(2:end) = curveSmo(2:end)>=curveSmo(1:end-1);
                right = false(numel(TW),1);
                right(end) = true;
                right(1:end-1) = curveSmo(1:end-1)>=curveSmo(2:end);
                maxima = left & right;
                maxima = find(maxima);
                maxima = maxima(curveSmo(maxima)>=maxV-3*s0);
                if(numel(maxima)==1)
                   tPeak = maxima + min(TW)-1;
                   break;
                end
                if(smo>=5)   % limitation on maximum smoothing
                   [~,id] = max(curveSmo(maxima));
                   tPeak = maxima(id) + min(TW)-1;
                    break;
                end
                smo = smo + 1;
            end
        else
            tPeak = TW;
        end

        if(tPeak == min(TW) || tPeak == max(TW)) 
            curveSmo = imgaussfilt(curve,2);
            [~,tPeak] = max(curveSmo);
            tPeak = tPeak + min(TW)  -1;
        end
        majorityEvt{i}.tPeak = tPeak;
        majorityEvt{i}.needUpdatePeak = false;
    end
    
    
    
    %% get neighbor relation
    Map = zeros([H,W,T],'uint16');
    for i=1:numel(evtLst)
        Map(evtLst{i}) = i;
    end
    
    dh = [-1 0 1 -1 1 -1 0 1];
    dw = [-1 -1 -1 0 0 1 1 1];
    neibLst = cell(N,1);
    exLstSpa = cell(N,1);
    delayDif = cell(N,1);
    evtCCLabel = zeros(N,1);
    labelsInActRegs = cell(numel(ccRegions),1);
    parfor i = 1:N
        pix = evtLst{i};
        [ih,iw,it] = ind2sub([H,W,T],pix);
        neib0 = [];
        % Connected edge
        for ii=1:numel(dh)
            ih1 = min(max(ih + dh(ii),1),H);
            iw1 = min(max(iw + dw(ii),1),W);
            vox1 = sub2ind([H,W,T],ih1,iw1,it);
            idxSel = setdiff(Map(vox1),[0,i]);
            neib0 = union(neib0,idxSel);
        end
        neibLst{i} = neib0;
    end
    
    %% get CC relation
    for iReg = 1:numel(ccRegions)
        labelsInActReg = setdiff(Map(ccRegions{iReg}),0);
        labelsInActRegs{iReg} = labelsInActReg;
        evtCCLabel(labelsInActReg) = iReg;
        if(size(labelsInActReg,2)==1)
            labelsInActReg = labelsInActReg';
        end
        %% update exLstSpa
        for curLabel = labelsInActReg
%             curLabel = labelsInActReg(i);
            ex0 = curLabel;
            curMajIhw =  majorityEvt{curLabel}.ihw;
            for j = 1:numel(labelsInActReg)
                nLabel = labelsInActReg(j);
                if(curLabel>=nLabel)
                    continue;
                end
                nMajIhw = majorityEvt{nLabel}.ihw;
                n0 = numel(intersect(curMajIhw,nMajIhw));
                n1 = numel(curMajIhw);
                n2 = numel(nMajIhw);
                if((n0/n1>overlap || n0/n2>overlap))
                    ex0 = [ex0;nLabel];
                end
            end
            exLstSpa{curLabel} = ex0;
        end
    end
    
    %% create delay map
    for i = 1:N
        delayDif{i} = containers.Map('KeyType','double','ValueType','double');
    end
    
    mergingInfo.neibLst = neibLst;
    mergingInfo.exLstSpa = exLstSpa;
    mergingInfo.delayDif = delayDif;
    mergingInfo.evtCCLabel = evtCCLabel;
    mergingInfo.labelsInActRegs = labelsInActRegs;
    mergingInfo.refineCheckList = true(N,1);
end