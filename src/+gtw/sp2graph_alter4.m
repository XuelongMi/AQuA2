function [ref,tst,refBase,s,t,idxGood,orgCurves,stds] = sp2graph_alter4(dFVec,spLst,major,spMap,WTW,spSz,seMap0,seSel)
% sp2graph convert super pixels to curves and graph nodes for GTW
% Input super pixels are not perfect:
% Empty or too small
% Bad corresponding curves
mIhw = major.ihw;
peakTW = major.TW;
tPeak = major.tPeak;
[H0,W0,T0] = size(spMap);
refBase = mean(dFVec(mIhw,:),1);

% make sure the maximum extension.
t00 = min(WTW);
seMap0 = reshape(seMap0,[],T0);
seMap0 = seMap0(mIhw,:);
for t = min(WTW)-1:-1:1
    labels = setdiff(seMap0(:,t),0);
    if(isempty(labels) || ismember(seSel,labels))
        t00 = t;
    else
        break; 
    end
end

ts = extendLeft(imgaussfilt(refBase,2),tPeak,min(peakTW),min(WTW),t00);
ext = 5;
validPart = false(T0,1);
validPart(ts:min(tPeak+ext,T0)) = true;
sz = cellfun(@numel,spLst);

%% make it non-decreasing
smo = 1;
while(1)
    curveSmo = imgaussfilt(refBase,smo);
    left = false(T0,1);
    left(1) = true;
    left(2:end) = curveSmo(2:end)>=curveSmo(1:end-1);
    right = false(T0,1);
    right(end) = true;
    right(1:end-1) = curveSmo(1:end-1)>=curveSmo(2:end);
    maxima = left & right & validPart;
    maxima = find(maxima);
    if(numel(maxima)<=1 || smo>=8)
       break;
    end
    smo = smo + 1;
end
ts = extendLeft(curveSmo,tPeak,min(peakTW),min(WTW),t00);

% smo = 2;curveSmo = imgaussfilt(refBase,smo);

refBase = curveSmo;
peakTW = ts:min(major.tPeak+ext,T0);
refBase(1:ts) = refBase(ts);
refBase(max(peakTW):end) = refBase(max(peakTW));
if(~isempty(maxima))
    [~,tPeak] = max(curveSmo(peakTW));
    tPeak = ts + tPeak - 1;
end
refBase = refBase - refBase(1);
refBase = refBase/refBase(tPeak);
refBase(tPeak:end) = 1;

%% non decreasing
refBase = monoIncreasing(refBase,tPeak,ts);

nSp = numel(spLst);
tst = zeros(nSp,T0);
ref = zeros(nSp,T0);
orgCurves = zeros(nSp,T0);
stds = zeros(nSp,1);

% ratio = sqrt(numel(mIhw)/spSz);
% if(ratio>1)
%     dF = reshape(dFVec,[H0,W0,T0]);
%     dF = imgaussfilt(dF,ratio/3);
%     dFVec = reshape(dF,[],T0);
% end

% smo = 2;
% valid = peakTW;  
mTW = major.TW;
mtPeak = major.tPeak;
estPeakShift = max(5,round(numel(major.TW)/3));
mtPeakTW = max(1,mtPeak-estPeakShift):min(mtPeak+estPeakShift,T0);
mTW = intersect(mTW,mtPeakTW);
valid = mTW;

for ii=1:numel(spLst)
    sp0 = spLst{ii};
    x0 = mean(dFVec(sp0,:),1);
    x0Smo = imgaussfilt(x0,2);
%     [vMax,tPeak] = max(x0Smo(valid));
%     tPeak = tPeak + min(valid) - 1;
    orgCurves(ii,:) = x0Smo;
    
    curTW = max(min(valid)-2,1):min(T0,max(valid)+2);
    left = false(numel(curTW),1);  left(2:end) = x0Smo(min(curTW)+1:max(curTW))>=x0Smo(min(curTW):max(curTW)-1);
    right = false(numel(curTW),1); right(1:end-1) = x0Smo(min(curTW):max(curTW)-1)>=x0Smo(min(curTW)+1:max(curTW));
    maxima = left & right; maxima = find(maxima); maxima = maxima + min(curTW) - 1;
%     [~,id] = max(x0Smo(maxima));
%     tPeak = maxima(id);
%     vMax = x0(tPeak);
%     if(isempty(maxima))
%         [~,id] = max(x0Smo(valid));
%         tPeak = valid(id);
% %         vMax = x0(tPeak);
%     end
    score = (x0Smo(maxima)').*exp(-(maxima-mtPeak).^2/2/estPeakShift^2/9);
    [~,id] = max(score);
    tPeak = maxima(id); if(isempty(maxima)) tPeak = mtPeak; end;
    
    %% cut at min
%     ts = extendLeft(x0,x0Smo,tPeak,min(peakTW),max(min(peakTW)-ext,1)); 
    
    ts1 = max(min(peakTW)-ext,1);
    [~,ts] = min(x0Smo(ts1:tPeak));
    ts = ts + ts1 - 1;
    
    stds(ii) =  mean((x0Smo(2:end)-x0Smo(1:end-1)).^2)/2;
    x0 = x0Smo;
    vMin = x0(ts);vMax = x0(tPeak);
%     x0(1:ts-1) = 0;
%     x0(tPeak+1:end) = 1;
%     x0(ts:tPeak) = (x0(ts:tPeak)-vMin)/(vMax-vMin);
%     tst(ii,:) = x0;
%     ref(ii,:) = refBase;
    
    x0(1:ts-1) = vMin;
    x0(tPeak+1:end) = vMax;
    x0 = x0 - vMin;
    x0 = monoIncreasing(x0,tPeak,ts);
    tst(ii,:) = x0;
    ref(ii,:) = refBase*(vMax-vMin);
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
function ts = extendLeft(x0Smo,tPeak,tMin,t00,t01)
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
function refBase = monoIncreasing(refBase,tPeak,ts)
    preV = refBase(tPeak);
    preT = tPeak;
    for t = tPeak:-1:ts
       if(refBase(t)<preV) 
           refBase(t:preT) = refBase(t) + (preV-refBase(t))/(preT-t)*[0:preT-t];
           preV = refBase(t);
           preT = t;
       end
    end
end