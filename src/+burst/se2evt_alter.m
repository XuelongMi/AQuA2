function [evtRecon,evtL,dlyMap,nEvt0,rgtx,svLabel] = se2evt_alter(...
        dF0,seMap0,seSel,ihw0,rgt,superVoxels,major0,opts)

gtwSmo = opts.gtwSmo; % 0.5

spSz = 25;  % preferred super pixel size  200
spT = 30;  % super pixel number scale (larger for more)

% GTW on super pixels
% group super pixels to events
[H,W,T] = size(dF0);
if numel(ihw0)>30
    [spLst,tDly,dlyMap,majoirtyMap,scl] = gtw.spgtw3(...
        dF0,seMap0,seSel,gtwSmo,spSz,superVoxels,major0);
    if(numel(superVoxels)>1)
        [svLabel] = burst.delayMap2evt(spLst,tDly,majoirtyMap,major0,spSz,opts);
    else
        svLabel = ones(numel(superVoxels),1);
    end
    
%     [spLst,tDly,neiLst,sv_spLabels,rgtSel,xFail,riseDur] = gtw.spgtw2(...
%         dF0,seMap0,seSel,gtwSmo,spSz,superVoxels,major0);
%     if xFail==0
        % smooth propagation first
%         [svLabel] = burst.riseMap2evtMuYu_Alter3(spLst,tDly,neiLst,sv_spLabels,maxRiseUnc,cDelay,[H,W],spSz,riseDur);
%     end
else
    scl = 1;
    rgtSel = 1:numel(rgt);
    spLst = {ihw0};
    svLabel = ones(numel(superVoxels),1);
%     dlyMap = evtMap;
    dlyMap = min(rgt)*ones(H,W);
    dF0Vec = reshape(dF0,[],numel(rgt));
end

% events
[H,W,T] = size(dF0);
evtL = zeros(size(dF0),'uint16');
for i = 1:numel(superVoxels)
    evtL(superVoxels{i}) = svLabel(i);
end
% rgtx = min(it0):max(it0);
[evtL,evtRecon] = gtw.evtRecon_alter(evtL,seMap0,dF0,seSel);
evtRecon = uint8(evtRecon*255);
evtInTemporalDimension = sum(reshape(evtL>0,[],T)>0,1);
t0 = find(evtInTemporalDimension>0,1);
t1 = find(evtInTemporalDimension>0,1,'last');
rgtx = rgt(t0:t1);
evtL = evtL(:,:,t0:t1);
evtRecon = evtRecon(:,:,t0:t1);
nEvt0 = max(svLabel);
dlyMap = dlyMap*scl + min(rgt)-1;
end

