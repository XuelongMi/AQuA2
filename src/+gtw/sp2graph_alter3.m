function [ref,tst,refBase,s,t,idxGood,orgCurves,stds] = sp2graph_alter3(dFVec,spLst,major,spMap,WTW)
% sp2graph convert super pixels to curves and graph nodes for GTW
% Input super pixels are not perfect:
% Empty or too small
% Bad corresponding curves
mIhw = major.ihw;
peakTW = major.TW;
tPeak = major.tPeak;
[H0,W0,T0] = size(spMap);
% spMap = reshape(spMap,[],T0);
refBase = mean(dFVec(mIhw,:),1);
ts = extendLeft(refBase,imgaussfilt(refBase,2),tPeak,min(peakTW),min(WTW),1);
peakTW = ts:max(peakTW);
refBase = imgaussfilt(refBase,2);
refBase(1:ts) = refBase(ts);
refBase(max(peakTW):end) = refBase(max(peakTW));

r1 = refBase;
% [~,ix] = max(r1);
bw = refBase*0;
bw(tPeak) = 1;
r2 = -imimposemin(-r1,bw);
r2(tPeak) = r1(tPeak);
r2(isinf(r2)) = nan;
refBase = r2;


[~,tPeak] = max(refBase);
refBase = refBase - refBase(1);
refBase = refBase/refBase(tPeak);
refBase(tPeak:end) = 1;

nSp = numel(spLst);
tst = zeros(nSp,T0);
ref = zeros(nSp,T0);
orgCurves = zeros(nSp,T0);
stds = zeros(nSp,1);
ext = 5;
for ii=1:numel(spLst)
    sp0 = spLst{ii};
    % the biggest component
%     BW = sum(spMap(sp0,:),1)>0;
%     cc = bwconncomp(BW);
%     cc = cc.PixelIdxList;
%     sz = cellfun(@numel,cc);
%     [~,id] = max(sz);
    valid = peakTW;%intersect(cc{id},TW);
    
%     t0 = valid(1);
    x0 = mean(dFVec(sp0,:),1);
    x0Smo = imgaussfilt(x0,2);
%     valid = TW;
    [vMax,tPeak] = max(x0Smo(valid));
    orgCurves(ii,:) = x0Smo;
%     tPeak = valid(tPeak);
    
    curTW = max(min(valid)-2,1):min(T0,max(valid)+2);
    left = false(numel(curTW),1);  left(2:end) = x0Smo(min(curTW)+1:max(curTW))>=x0Smo(min(curTW):max(curTW)-1);
    right = false(numel(curTW),1); right(1:end-1) = x0Smo(min(curTW):max(curTW)-1)>=x0Smo(min(curTW)+1:max(curTW));
    maxima = left & right; maxima = find(maxima); maxima = maxima + min(curTW) - 1;
    [~,id] = max(x0Smo(maxima));
    tPeak = maxima(id);
%     vMax = x0(tPeak);
    if(isempty(maxima))
        [~,id] = max(x0Smo(valid));
        tPeak = valid(id);
%         vMax = x0(tPeak);
    end
    
    %% cut at min
%     ts = extendLeft(x0,x0Smo,tPeak,min(peakTW),max(min(peakTW)-ext,1)); 
    
    ts1 = max(min(peakTW)-ext,1);
    [~,ts] = min(x0Smo(ts1:tPeak));
    ts = ts + ts1 - 1;
    
    stds(ii) =  mean((x0Smo(2:end)-x0Smo(1:end-1)).^2)/2;
    x0 = x0Smo;
    vMin = x0(ts);vMax = x0(tPeak);
    x0(1:ts-1) = 0;
    x0(tPeak+1:end) = 1;
    x0(ts:tPeak) = (x0(ts:tPeak)-vMin)/(vMax-vMin);
    tst(ii,:) = x0;
    ref(ii,:) = refBase;
    %     end
end

idxGood = var(tst,0,2)>1e-10;
% spLst = spLst(idxGood);
% tst = tst(idxGood,:);
% ref = ref(idxGood,:);

% graph, at most one pair between two nodes
s = nan(nSp,1);
t = nan(nSp,1);
nPair = 0;
dh = [-1 0 1 -1 0 1 -1 0 1];
dw = [-1 -1 -1 0 0 0 1 1 1];
% tEnd = min(T0,tEnd+5);
% tst = tst(:,1:tEnd);
% ref = ref(:,1:tEnd);
% refBase = refBase(:,1:tEnd);
spMap1 = zeros(H0,W0);
for ii=1:numel(spLst)
    spMap1(spLst{ii}) = ii;
end
nSp = numel(spLst);
for ii=1:nSp
    sp0 = spLst{ii};
    [ih0,iw0] = ind2sub([H0,W0],sp0);
    neib0 = [];
    for jj=1:numel(dh)  % find neighbors in eight directions
        ih = max(min(H0,ih0+dh(jj)),1);
        iw = max(min(W0,iw0+dw(jj)),1);
        ihw = sub2ind([H0,W0],ih,iw);
        neib0 = union(neib0,spMap1(ihw));
    end
    neib0 = neib0(neib0>ii);    
    ids = nPair+[1:numel(neib0)];
    nPair = nPair + numel(neib0);
    s(ids) = ii;
    t(ids) = neib0;
end
s = s(~isnan(s));
t = t(~isnan(t));

end
function [curve,tPeak] = balanceCurve(curve)
    [vP,tPeak] = max(curve);
    [vM,tMin] = min(curve(1:tPeak));
    curve(1:tMin) = vM;
    curve(tPeak:end) = vP;
    curve = curve - vM;
    curve = curve/max(curve);
end
function ts = extendLeft(x0,x0Smo,tPeak,tMin,t00,t01)
%     x0 = x0Smo;
% 	s0 = mean((x0Smo(2:end)-x0Smo(1:end-1)).^2)/2;
    
    [vMin,~] = min(x0Smo(t00:tPeak));
    thr = 0.3*(x0Smo(tPeak)-vMin);
%     tMin = tPeak; vMin = x0Smo(tPeak);
    ts = tMin;
    vMin = x0Smo(tMin);
    
    %% extend backward
    for t = tMin:-1:t01
        if(x0Smo(t)-vMin>=thr)
            break;
        else
           if(x0Smo(t)<vMin) 
              vMin = x0Smo(t) ;
              ts = t;
           end
        end
    end
end