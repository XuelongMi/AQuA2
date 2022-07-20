function [arLst] = acDetect(dF,opts,evtSpatialMask,ff)
    thrs = opts.maxdF1:(opts.thrARScl-opts.maxdF1)/20:opts.thrARScl;
    [H,W,T] = size(dF);
    evtSpatialMask = evtSpatialMask(:);
    dF = reshape(dF,[],T);
    dF(~evtSpatialMask,:) = 0;
    dF = reshape(dF,[H,W,T]);
    % valid region
    activeMap = false(H,W,T);
    tic;
    for k = 1:numel(thrs)
        thr = thrs(k);
        if(exist('ff','var'))
            waitbar(0.5 + 0.5*k/numel(thrs),ff,'Detecting active regions...');
        end
        selectMap = (dF>thr);
        curRegions = bwconncomp(selectMap);
        curRegions = curRegions.PixelIdxList;
        valid = false(numel(curRegions),1);
        parfor i = 1:numel(curRegions)
            pix = curRegions{i};
            [ih,iw,it] = ind2sub([H,W,T],pix);
            rgh = min(ih):max(ih); H0 = numel(rgh); ih = ih - min(ih) + 1;
            rgw = min(iw):max(iw); W0 = numel(rgw); iw = iw - min(iw) + 1;
            rgt = min(it):max(it); T0 = numel(rgt); it = it - min(it) + 1;
            pix0 = sub2ind([H0,W0,T0],ih,iw,it);
            curMap = zeros(H0,W0,T0,'single');
            curMap(pix0) = 1;
            curMap = (sum(curMap,3)>opts.compress*T0);
            ihw = find(curMap);
            curSz = numel(ihw);
            
            % size, duration limitation
            if (curSz>opts.maxSize || curSz<opts.minSize || T0<opts.minDur)
               continue; 
            end

            erodeMap = imerode(curMap,strel('disk',1));
            boundary = curSz - sum(erodeMap(:));
            circularity = 4*pi*curSz/(boundary^2);
            
            % circularity limitation
            if(circularity>opts.circularityThr)
                valid(i) = true;
            end
        end
        curRegions = curRegions(valid);
        for i = 1:numel(curRegions)
            activeMap(curRegions{i}) = true;
        end
    end
    toc;
    arLst = bwconncomp(activeMap);
    arLst = arLst.PixelIdxList;
end

