function mergingInfo = gapIdentify(evtLst,dF,dFOrg,mergingInfo,opts)
    [H,W,T] = size(dF);
%     eig3 = cal_pc_3D(double(dF),5);
%     gapThr = 0.5;
    %% setting
    round = 10;
    gapSigThr = 50;
    [dw,dh,dt] = se.dirGenerate(26);
    
    %% normalize dFOrg
    xx = (dFOrg(:,:,2:end) - dFOrg(:,:,1:end-1)).^2;
    stdMap = sqrt(median(xx,3)/0.9133);
    stdMapGau = double(imgaussfilt(stdMap)) + 1e-6;
    dFOrg = dFOrg./repmat(stdMapGau,1,1,T) + opts.xBias;
    
    %%
    evtLstGrow = cell(numel(evtLst),1);
    for i = 1:numel(evtLst)
        [ih0,iw0,it0] = ind2sub([H,W,T],evtLst{i});
        curGrow = [];
        for k = 1:numel(dw)
           ih = max(1,min(H,ih0+dh(k)));
           iw = max(1,min(W,iw0+dw(k)));
           it = max(1,min(T,it0+dt(k)));
           newAdd = sub2ind([H,W,T],ih,iw,it);
           curGrow = [curGrow;newAdd];
        end
        evtLstGrow{i} = unique(curGrow);        
    end
    
    for i = 1:numel(mergingInfo.neibLst)
        curLabel = i;
        curGrow = evtLstGrow{curLabel};
        neiLst = mergingInfo.neibLst{i};
        neiLst = neiLst(neiLst>curLabel);
        checkGap = true(numel(neiLst),1);
        maxCurV = max(dF(evtLst{curLabel}));
        for j = 1:numel(neiLst)
            nLabel = neiLst(j);
            nGrow = evtLstGrow{nLabel};
            inter = intersect(curGrow,nGrow);
            maxnV = max(dF(evtLst{nLabel}));
%             if(mean(eig3(inter))>gapThr)
%                 checkGap(j) = false;
%             end
            
            if(3*mean(dF(inter))<min(maxCurV,maxnV))
                checkGap(j) = false;
            end

%             %% grow boundary
%             growRound = cell(round*2,1);
%             curGrow = inter;
%             growRound{1} = inter;
%             Range = union(evtLst{curLabel},evtLst{nLabel});
%             for k = 2:round*2
%                curGrow = growOneRound(curGrow,dw,dh,dt,H,W,T,Range);
%                growRound{k} = curGrow;
%             end
%             %% order statistic + T-test
%             for k = 1:round
%                 bg = growRound{k};
%                 fg = setdiff(growRound{2*k},bg);
%                 if(isempty(fg)||isempty(bg))
%                     continue;
%                 end
%                 fg = dFOrg(fg); bg = dFOrg(bg);
%                 L = mean(fg)-mean(bg);
%                 [mu, sigma] = ksegments_orderstatistics_fin(fg, bg);
%                 zscore = (L-mu)/sigma;
%                 T_zscore = (mean(fg)-mean(bg))/sqrt(1/numel(fg)+1/numel(bg));
%                 if(zscore>gapSigThr && T_zscore>gapSigThr)
%                     checkGap(j) = false;
%                     break;
%                 end
%             end
        end
        neiLst = neiLst(checkGap);
        mergingInfo.neibLst{i} = neiLst;
    end

end
function pix = growOneRound(pix,dw,dh,dt,H,W,T,Range)

    [ih0,iw0,it0] = ind2sub([H,W,T],pix);
    curGrow = [];
    for k = 1:numel(dw)
       ih = max(1,min(H,ih0+dh(k)));
       iw = max(1,min(W,iw0+dw(k)));
       it = max(1,min(T,it0+dt(k)));
       newAdd = sub2ind([H,W,T],ih,iw,it);
       curGrow = [curGrow;newAdd];
    end
    pix = intersect(Range,curGrow);
end