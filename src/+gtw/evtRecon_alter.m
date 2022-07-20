function [evtL,evtRecon] = evtRecon_alter(evtL,seMap0,dF0,seSel)

evtLst = label2idx(evtL);
[H,W,T] = size(seMap0);
evtL = reshape(evtL,[],T);
evtRecon = zeros(H*W,T);
evtInTemporalDimension = sum(evtL>0,1);
t0 = find(evtInTemporalDimension>0,1);
t1 = find(evtInTemporalDimension>0,1,'last');

spSz = 25;
nSpEstimate = max(1,round(H*W/spSz));
dFAvg = max(dF0(:,:,t0:t1),[],3);
L = superpixels(dFAvg,nSpEstimate,'Compactness',20);
spLst = label2idx(L);
dF0 = imgaussfilt(dF0,2);
dF0Vec = reshape(dF0,[],T);
seMap0 = reshape(seMap0,[],T);

whetherExtend = 0;

for ii=1:numel(evtLst)
    pix = evtLst{ii};
    [ih,iw,it] = ind2sub([H,W,T],pix);
    ihw = unique(sub2ind([H,W],ih,iw));
    spLabels = setdiff(L(ihw),0);
    evtL0 = evtL==ii;
    for k = 1:numel(spLabels)
        curSp = spLst{(spLabels(k))};
        curSp = intersect(curSp,ihw);
        x0 = mean(dF0Vec(curSp,:),1);
        s0 = mean((x0(2:end)-x0(1:end-1)).^2)/2;
        curTemporal = sum(evtL0(curSp,:),1);
        t0 = find(curTemporal>0,1);
        t1 = find(curTemporal>0,1,'last');
        maxV = max(x0(t0:t1));
        minV = x0(t0);ts = t0;
        % left
        for t = t0:-1:1
            curse = seMap0(curSp,t); cure = evtL(curSp,t);
             curse = (curse>0 & curse~=seSel) | (cure>0 & cure~=ii);
            if(sum(curse(:))==numel(curSp))
                break;
            end
            if(x0(t)<minV)
                minV = x0(t);
                ts = t;
            else
                if(x0(t)<0.1*maxV || x0(t)-minV>=3*s0)
                    break;
                end
            end
        end
        % right
        minV = x0(t1); te = t1;
        for t = t1:T
            curse = seMap0(curSp,t); cure = evtL(curSp,t);
            curse = (curse>0 & curse~=seSel) | (cure>0 & cure~=ii);
            if(sum(curse(:))==numel(curSp))
                break;
            end
            if(x0(t)<minV)
                minV = x0(t);
                te = t;
            else
                if(x0(t)<0.1*maxV || x0(t)-minV>=3*s0)
                    break;
                end
            end
        end
        
        % update evtL and evtRecon
        if(whetherExtend == 1)
            avaliable = seMap0(curSp,ts:te)==0 |evtL(curSp,ts:te)==0 | evtL(curSp,ts:te)==ii;
        else
            avaliable = evtL(curSp,ts:te)==ii;
        end
        recon = repmat(x0/maxV,numel(curSp),1);
        recon = recon(:,ts:te);
        curEvtL = evtL(curSp,ts:te);
        curEvtL(avaliable) = ii;
        evtL(curSp,ts:te) = curEvtL;
        curRecon = evtRecon(curSp,ts:te);
        curRecon(avaliable) = recon(avaliable);
        evtRecon(curSp,ts:te) = curRecon;
    end
end
evtL = reshape(evtL,H,W,T);
evtRecon = reshape(evtRecon,H,W,T);
end
