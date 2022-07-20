function [subEvtLst,seLabel,majorInfo] = preSplitSub(subEvtLst,seLabel,majorInfo,sz)
    H = sz(1);W = sz(2);T = sz(3);
    nEvt = numel(subEvtLst);
    for i = 1:numel(subEvtLst)
        pix = subEvtLst{i};
        [ih,iw,it] = ind2sub([H,W,T],pix);
        rgh = min(ih):max(ih); H0 = numel(rgh);
        rgw = min(iw):max(iw); W0 = numel(rgw);
        rgt = min(it):max(it); T0 = numel(rgt);
        Map = false(H0,W0,T0);
        ih = ih - min(rgh) + 1;
        iw = iw - min(rgw) + 1;
        it = it - min(rgt) + 1;
        pix0 = sub2ind([H0,W0,T0],ih,iw,it);
        Map(pix0) = true;
        cc = bwconncomp(Map);
        if(cc.NumObjects>1)
            cc = cc.PixelIdxList;
            for j = 1:numel(cc)
                if(j>1)
                    nEvt = nEvt + 1;
                    curLabel = nEvt;
                else
                    curLabel = i;
                end
                pix0 = cc{j};
                [ih,iw,it] = ind2sub([H0,W0,T0],pix0);
                ih = ih + min(rgh) - 1;
                iw = iw + min(rgw) - 1;
                it = it + min(rgt) - 1;
                pix = sub2ind([H,W,T],ih,iw,it);
                subEvtLst{curLabel} = pix;
                seLabel(curLabel) = seLabel(i);
                majorInfo{curLabel} = majorInfo{i};
            end
        end
    end
end