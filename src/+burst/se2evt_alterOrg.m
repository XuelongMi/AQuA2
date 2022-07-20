function [evtRecon,evtL,dlyMap,nEvt0,rgtx] = se2evt_alterOrg(...
        dF0,seMap0,seSel,ihw0,rgt,superVoxels,major0,opts)

gtwSmo = opts.gtwSmo; % 0.5
maxStp = opts.maxStp; % 11
maxRiseUnc = opts.cRise;  % 1
cDelay = opts.cDelay;  % 5

spSz = 100;  % preferred super pixel size  200
spT = 30;  % super pixel number scale (larger for more)

% GTW on super pixels
% group super pixels to events
[H,W,T] = size(dF0);
if numel(ihw0)>30    
    [spLst,tDly,neiLst,sv_spLabels,rgtSel,xFail,riseDur] = gtw.spgtw2(...
        dF0,seMap0,seSel,gtwSmo,spSz,superVoxels,major0);
    if xFail==0
%         smooth propagation first
        [svLabel] = burst.riseMap2evtMuYu_Alter3(spLst,tDly,neiLst,sv_spLabels,maxRiseUnc,cDelay,[H,W],spSz,riseDur);
    end
end
    
if numel(ihw0)<=30 || xFail==1
    rgtSel = 1:numel(rgt);
    spLst = {ihw0};
    svLabel = ones(numel(superVoxels),1);
%     dlyMap = evtMap;
    dlyMap = [];
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
dlyMap = [];
end

