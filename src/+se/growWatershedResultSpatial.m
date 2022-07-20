function [evtLst,ccRegions] = growWatershedResultSpatial(evtLst,majorityEvt,dF,opts)
    [H,W,T] = size(dF);
    
    dw = [-1,-1,-1,0,0,1,1,1];
    dh = [-1,0,1,-1,1,-1,0,1];
    opts.slackThr = 1;
    Map = zeros(H,W,T,'uint16');
    Map(dF<opts.thrARScl - opts.slackThr) = 65535;
    dF = reshape(dF,[],T);
%     Map(~activeMap) = 65535;
    evt2D = cell(numel(evtLst),1);
    evt2D_cantAdd = cell(numel(evtLst),1);
    new_add_Pix = cell(numel(evtLst),1);
    for i = 1:numel(evtLst)
        Map(evtLst{i}) = i;
        pix = evtLst{i};
        [ih0,iw0,it0] = ind2sub([H,W,T],pix);
        ihw = unique(sub2ind([H,W],ih0,iw0));
        evt2D{i} = ihw;
        new_add_Pix{i} = ihw;
    end
    Map = reshape(Map,[],T);
    
    for k = 1:40
%         tic;
        disp(['Grow Round: ',num2str(k)]);
        for i = 1:numel(evtLst)
            ihw = new_add_Pix{i};
            [ih0,iw0] = ind2sub([H,W],ihw);
            new_add = [];
            for ii = 1:numel(dw)
               ih = max(1,min(H,ih0+dh(ii)));
               iw = max(1,min(W,iw0+dw(ii)));
               new_add = union(new_add,unique(sub2ind([H,W],ih,iw)));
            end
            new_add = setdiff(new_add,evt2D{i} );
            new_add = setdiff(new_add,evt2D_cantAdd{i});
            
            if(isempty(new_add))
                continue;
            end
            
            TW = majorityEvt{i}.TW;
            mIhw = majorityEvt{i}.ihw;
            curve = mean(dF(mIhw,:),1);
            if(numel(TW)<3)
               [maxV,tPeak] = max(curve(TW));
               tPeak = tPeak + min(TW) - 1;
               t0 = tPeak;
               t1 = tPeak;
               lMin = maxV; rMin = maxV;
               for t = tPeak-1:-1:1
                   if(curve(t)<lMin)
                       lMin = curve(t);
                   else
                       t0 = t+1;
                      break; 
                   end
               end
               
               for t = tPeak+1:T
                   if(curve(t)<rMin)
                       rMin = curve(t);
                   else
                       t1 = t-1;
                      break; 
                   end
               end
               TW = t0:t1;
            end
            [~,tPeak] = max(curve(TW));
            curve = curve(TW);
            curve0 = dF(new_add,TW);
            r = corr(curve',curve0');
            zscore = r./sqrt(1-r.^2)*sqrt(numel(TW)-2);
            pVal = tcdf(zscore,numel(TW)-2,'upper');
            if(size(new_add,1)==1)
                new_add = new_add';
            end
            select = pVal<1e-3;
            evt2D_cantAdd{i} = [evt2D_cantAdd{i};new_add(~select);];
            new_add = new_add(select);
            if(isempty(new_add))
                new_add_Pix{i} = [];
                continue;
            end
            
            curve0 = curve0(select,:);
            valid_add = false(numel(new_add),1);
            for j = 1:numel(new_add)
                lblCur = Map(new_add(j),TW);
                t0 = find(lblCur(1:tPeak)>0,1,'last');
                t0 = t0 + 1;
                if(t0>tPeak)
                    continue;
                end
                t1 = find(lblCur(tPeak:end)>0,1);
                t1 = t1 + tPeak - 1 - 1;
                if(t1<tPeak)
                    continue;
                end
                if isempty(t0)
                    t0 = 1;
                end
                if isempty(t1)
                    t1 = numel(lblCur);
                    
                end
                [leftMin,t_l] = min(curve0(j,t0:tPeak));
                t0 = t0 + t_l - 1;
                [rightMin,t_r] = min(curve0(j,tPeak:t1));
                t1 = tPeak + t_r - 1;
                maxV = curve0(j,tPeak);
                lblCur(t0:t1) = i;
                isGood = sum(lblCur==i)>0;
                if(isGood)
                    Map(new_add(j),TW) = lblCur;
                    valid_add(j) = true;
                end
            end
            evt2D{i} = [evt2D{i};new_add(valid_add)];
            evt2D_cantAdd{i} = [evt2D_cantAdd{i};new_add(~valid_add)];
            new_add_Pix{i} = new_add(valid_add);
        end
%         toc;
    end
    Map(Map==65535) = 0;
    Map = reshape(Map,[H,W,T]);
    evtLst = label2idx(Map);
    Map = zeros(H,W,T);
    for i = 1:numel(evtLst)
        pix = evtLst{i};
        [ih,iw,it] = ind2sub([H,W,T],pix);
        rgH = min(ih):max(ih);
        rgW = min(iw):max(iw);
        rgT = min(it):max(it);
        H0 = numel(rgH);
        W0 = numel(rgW);
        T0 = numel(rgT);
        ih = ih - min(ih) + 1;
        iw = iw - min(iw) + 1;
        it = it - min(it) + 1;
        pix = sub2ind([H0,W0,T0],ih,iw,it);
        Map0 = false(numel(rgH),numel(rgW),numel(rgT));
        Map0(pix) = true;
        curR = bwconncomp(Map0);
        curR = curR.PixelIdxList;
        sz = cellfun(@numel,curR);
        [~,id] = max(sz);
        pix = curR{id};
        [ih,iw,it] = ind2sub([H0,W0,T0],pix);
        ih = ih + min(rgH) - 1;
        iw = iw + min(rgW) - 1;
        it = it + min(rgT) - 1;
        pix = sub2ind([H,W,T],ih,iw,it);
        evtLst{i} = pix;
        Map(pix) = i;
    end
    
    ccRegions = bwconncomp(Map>0);
    ccRegions = ccRegions.PixelIdxList;


end