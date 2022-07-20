function [evtLst,sdLst,curRegions] = markerControlledSplitting_Ac(Map,activeMap,dF,dFOrg,opts)
    [H,W,T] = size(dF);
    sdLst = label2idx(Map);
    curThr = opts.thrARScl;
    selectMap = (dF>curThr) & activeMap;
    scoreMap = -imgaussfilt(dF,opts.spaSmo);% spatial smoothing for weakening gap in spatial
    clear dF;
    curRegions = bwconncomp(selectMap);
    curRegions = curRegions.PixelIdxList;
    dFVec = reshape(dFOrg,[],T);
    clear dFOrg;
    % calculate scoreMap
    SE = strel('cube',8);
    Ttestopen = true;
    nEvt = max(Map(:));
    gap = 2*opts.smoT+1;
    seedsInRegion = cell(numel(curRegions),1);
    % check whether the whole active region is significant or not
    for i = 1:numel(curRegions)
         pix = curRegions{i};
         [ih,iw,it] = ind2sub([H,W,T],pix);
         ihw = unique(sub2ind([H,W],ih,iw));
         dur = max(it)-min(it)+1;
         curRes = [];
         labels = setdiff(Map(pix),0);
         
         if(isempty(labels))    % no seed, check
            curLoc = [];
            curve = mean(dFVec(ihw,:),1);
            s0 = sqrt(median((curve(gap+1:end)-curve(1:end-gap)).^2)/0.9133);
            curve = curve/s0;
            t00 = min(it);
            t11 = max(it);
            stopCheck = false;
            if(dur<opts.minDur || numel(ihw)<opts.minSize)
                stopCheck = true;
             end 
            thrs = max(curve(t00:t11)):-opts.ratio:floor(min(curve(t00:t11)));
            for k = 1:numel(thrs)
                if(stopCheck)
                   break; 
                end
               thr = thrs(k);
               curTWs = bwconncomp(curve(t00:t11)>thr);
               curTWs = curTWs.PixelIdxList;
               sz = cellfun(@numel,curTWs);
               curTWs = curTWs(sz>=opts.minDur);
               for j = 1:numel(curTWs)
                   TW = curTWs{j} + t00 - 1;
                   t0 = TW(1);
                   t1 = TW(end);
                   fg = curve(TW);
                   [bg1,bg2] = se.findNeighbor(curve,t0,t1,T);
                   [Tscore1,Tscore2] = se.calTScore(fg,bg1,bg2);
                   if(Ttestopen && (Tscore1<=opts.sigThr || Tscore2<=opts.sigThr)) % T-test check
                      continue; 
                   end
                   [score1,score2] = se.calOrderScore(fg,bg1,bg2);
                   if(score1>opts.sigThr && score2>opts.sigThr)    % Update
                       nEvt = nEvt+1;
                       Map(pix) = nEvt;
                       sdLst{nEvt} = pix;
                       stopCheck = true;
                       labels = nEvt;
                   end
               end   
            end
         elseif(numel(labels) == 1)
            Map(pix) = labels(1);
         end
         seedsInRegion{i} = labels;
    end
            
    % watershed
     splitRegion = cell(numel(curRegions),1);
     for i = 1:numel(curRegions)
        labels = seedsInRegion{i};
        if(numel(labels)>1)
            pix = curRegions{i};
            [ih,iw,it] = ind2sub([H,W,T],pix);
            % Multiple seeds, need to split 
            rgh = min(ih):max(ih); H0 = numel(rgh); ih = ih - min(ih) + 1;
            rgw = min(iw):max(iw); W0 = numel(rgw); iw = iw - min(iw) + 1;
            rgt = min(it):max(it); T0 = numel(rgt); it = it - min(it) + 1;
            pix0 = sub2ind([H0,W0,T0],ih,iw,it);
            % Map seed regions
            Map0 = zeros(H0,W0,T0,'uint16');
            Map0(pix0) = Map(pix);
            % watershed input
            scoreMap0 = scoreMap(rgh,rgw,rgt);
            % deal with background
            BW = true(size(Map0));
            BW(pix0) = false;
            BW2 = imerode(BW,SE);
            scoreMap0(BW) =  max(scoreMap(pix));
            scoreMap0(BW2) = 0;
            scoreMap1 = imimposemin(scoreMap0,Map0>0|BW2,26);
            % marker-controlled splitting
            MapOut = watershed(scoreMap1,26);
            % update
            MapOut(BW) = 0;
            waterLst = label2idx(MapOut);
            curLoc = cell(numel(labels),1);
            for ii = 1:numel(labels)
               target = sdLst{labels(ii)};
               [ih0,iw0,it0] = ind2sub([H,W,T],target(1));
               ih0 = ih0 - min(rgh) + 1;
               iw0 = iw0 - min(rgw) + 1;
               it0 = it0 - min(rgt) + 1;

               waterLabels = MapOut(ih0,iw0,it0);
               curPix = [waterLst{waterLabels}];
               [ch,cw,ct] = ind2sub([H0,W0,T0],curPix);
               ch = ch + min(rgh) - 1;
               cw = cw + min(rgw) - 1;
               ct = ct + min(rgt) - 1;
               curLoc{ii} = sub2ind([H,W,T],ch,cw,ct);
            end
            splitRegion{i} = curLoc;
        end
     end
     
    % update Map, grow one circle to remove gap
    [x_dir,y_dir,z_dir] = se.dirGenerate(26); 
    validRegion = false(numel(curRegions),1);
    for i = 1:numel(curRegions)
        labels = seedsInRegion{i};
        if(~isempty(labels))
            validRegion(i) = true;
        end
        if(numel(labels)>1)
           for ii = 1:numel(labels)
               pix = splitRegion{i}{ii};
               Map(pix) = labels(ii);
           end
           pix = curRegions{i};
           pix = pix(Map(pix)==0);
           [ih0,iw0,it0] = ind2sub([H,W,T],pix);
           for k = 1:numel(x_dir)
                ih = max(1,min(H,ih0+x_dir(k)));
                iw = max(1,min(W,iw0+y_dir(k)));
                it = max(1,min(T,it0+z_dir(k)));
                pixCur = sub2ind([H,W,T],ih,iw,it);
                select = Map(pixCur)>0;
                Map(pix(select)) = Map(pixCur(select));
                ih0 = ih0(~select);
                iw0 = iw0(~select);
                it0 = it0(~select);
                pix = pix(~select);
           end
        end
    end
    evtLst = label2idx(Map);
    curRegions = curRegions(validRegion);
end
