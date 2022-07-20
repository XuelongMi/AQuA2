function [riseLst,datR,evtLst,seLst] = se2evtTop(dF,seLst,svLst,seLabel,majorInfo,opts,ff)
    % evtTop super voxels to super events and optionally, to events
    load('random_Seed.mat');
    rng(s);
    gaptxx = opts.gapExt;
    [H,W,T] = size(dF);
    seMap = zeros(H,W,T,'uint16');
    for i = 1:numel(seLst)
       seMap(seLst{i})  = i;
    end
    
    % super event to events
    fprintf('Detecting events ...\n')
    riseLst = cell(0);
    datR = zeros(H,W,T,'uint8');
    datL = zeros(H,W,T,'uint16');
    nEvt = 0;
    for nn=  1:numel(seLst)
        se0 = seLst{nn};
        if isempty(se0)
            continue
        end
        fprintf('SE %d \n',nn)
        if exist('ff','var')&& ~isempty(ff)
            waitbar(0.1+nn/numel(seLst)*0.9,ff);
        end
        
        % super event pixel transform
        [ih0,iw0,it0] = ind2sub([H,W,T],se0);
        rgh = min(ih0):max(ih0); rgw = min(iw0):max(iw0);
        ihw0 = unique(sub2ind([numel(rgh),numel(rgw)],ih0-min(rgh)+1,iw0-min(rgw)+1));
        gapt = min(max(it0)-min(it0),gaptxx);rgt = max(min(it0)-gapt,1):min(max(it0)+gapt,T);
        
        % sub event pixel transform
        svLabels = find(seLabel==nn);
        superVoxels = cell(numel(svLabels),1);
        major0 = majorInfo(svLabels);
        for k = 1:numel(svLabels)
           pix =  svLst{svLabels(k)};
           [ih,iw,it] = ind2sub([H,W,T],pix);
           ihw = unique(sub2ind([H,W],ih,iw));
           ih = ih - min(rgh) + 1;
           iw = iw - min(rgw) + 1;
           it = it - min(rgt) + 1;
           superVoxels{k} = sub2ind([numel(rgh),numel(rgw),numel(rgt)],ih,iw,it);
           TW0 = major0{k}.TW;
           TW0 = max(min(rgt),min(TW0)):min(max(rgt),max(TW0));           
           TW0 = TW0 - min(rgt) + 1;

           major0{k}.TW = TW0;
           major0{k}.tPeak = major0{k}.tPeak - min(rgt) + 1;
           mIhw = major0{k}.ihw;
           if(numel(mIhw)<opts.minSize)
               mIhw = [];
           end
           [mIh,mIw] = ind2sub([H,W],mIhw);
           mIh = mIh - min(rgh) + 1;
           mIw = mIw - min(rgw) + 1;
           major0{k}.ihw = sub2ind([numel(rgh),numel(rgw)],mIh,mIw);
        end
        
        dF0 = dF(rgh,rgw,rgt);
        seMap0 = seMap(rgh,rgw,rgt);
         [evtRecon,evtL,dlyMap,nEvt0,rgtx,svLabel] = burst.se2evt_alter(...
            dF0,seMap0,nn,ihw0,rgt,superVoxels,major0,opts);
%         [evtRecon,evtL,dlyMap,nEvt0,rgtx] = burst.se2evt_alterOrg(...
%             dF0,seMap0,nn,ihw0,rgt,superVoxels,major0,opts);
        
        seMap00 = seMap(rgh,rgw,rgtx);
        evtL(seMap00~=nn & seMap00>0) = 0;
        evtL(evtL>0) = evtL(evtL>0)+nEvt;
        dLNow = datL(rgh,rgw,rgtx);
        dRNow = datR(rgh,rgw,rgtx);
        ixOld = dLNow>0;
        evtL(ixOld) = dLNow(ixOld);
        datR(rgh,rgw,rgtx) = max(dRNow,evtRecon);  % combine events
        datL(rgh,rgw,rgtx) = evtL;
        riseLst = burst.addToRisingMap(riseLst,superVoxels,svLabel,dlyMap,nEvt,rgh,rgw,rgt);
        nEvt = nEvt + nEvt0;
    end
    
    evtLst = label2idx(datL);    
end



