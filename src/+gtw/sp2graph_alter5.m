function [ref,tst,refBase,s,t,idxGood,orgCurves,stds] = sp2graph_alter5(rDF,df0ip,m0Msk,rgt00,vMap0,spLst,seedIn,gapSeedHW,majorInfo,sv_spLabels)
% sp2graph convert super pixels to curves and graph nodes for GTW
% Input super pixels are not perfect:
% Empty or too small
% Bad corresponding curves

nSp = numel(spLst);
sp_svLabels = zeros(nSp,1);
for i = 1:numel(sv_spLabels)
    sp_svLabels(sv_spLabels{i}) = i;
end
df0ip = df0ip(:,:,rgt00);
m0Msk = m0Msk(:,:,rgt00);
[H0,W0,T0] = size(df0ip);
[ih,iw] = ind2sub([H0,W0],seedIn);
rgh = max(ih-gapSeedHW,1):min(ih+gapSeedHW,H0);
rgw = max(iw-gapSeedHW,1):min(iw+gapSeedHW,W0);
df00Vec = reshape(rDF(rgh,rgw,:),[],T0);
vm00 = vMap0(rgh,rgw);
df00Vec = df00Vec(vm00>0,:);
refBase = nanmean(df00Vec,1);
% refBase = imgaussfilt(refBase,1);

[~,tPeak] = max(imgaussfilt(refBase,2));
ext = 5; ts = 1;
validPart = false(T0,1);
validPart(ts:min(tPeak+ext,T0)) = true;
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
    if(numel(maxima)<=1 || smo>=5)
       break;
    end
    smo = smo + 1;
end
refBase = curveSmo;
%%
peakTW = ts:min(tPeak,T0);
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
tEnd = 0;
preV = 1;
preT = tPeak;
for t = tPeak:-1:ts
   if(refBase(t)<preV) 
       refBase(t:preT) = refBase(t) + (preV-refBase(t))/(preT-t)*[0:preT-t];
       preV = refBase(t);
       preT = t;
   end
end

df0Vec = reshape(df0ip,[],T0);
m0Msk = reshape(m0Msk,[],T0);

nSp = numel(spLst);
tst = zeros(nSp,T0);
ref = zeros(nSp,T0);
orgCurves = zeros(nSp,T0);

estPeakShift = 10;
stds = zeros(numel(spLst),1);
for ii=1:numel(spLst)
    sp0 = spLst{ii};
    mTW = majorInfo{sp_svLabels(ii)}.TW - min(rgt00)+1;
    mtPeak = majorInfo{sp_svLabels(ii)}.tPeak - min(rgt00)+1;
    mtPeakTW = max(1,mtPeak-estPeakShift):min(mtPeak+estPeakShift,T0);
    mTW = intersect(mTW,mtPeakTW);
    valid = intersect(find(sum(m0Msk(sp0,:))>0),mTW);
    if(isempty(valid)) valid = mTW; end;
    t0 = valid(1);
    x0Org = mean(df0Vec(sp0,:),1);
    x0 = imgaussfilt(x0Org,2);
    [vMax,tPeak] = max(x0(valid));
    orgCurves(ii,:) = x0;
%     tPeak = valid(tPeak);
    
    curTW = max(min(valid)-2,1):min(T0,max(valid)+2);
    left = false(numel(curTW),1);  left(2:end) = x0(min(curTW)+1:max(curTW))>=x0(min(curTW):max(curTW)-1);
    right = false(numel(curTW),1); right(1:end-1) = x0(min(curTW):max(curTW)-1)>=x0(min(curTW)+1:max(curTW));
    maxima = left & right; maxima = find(maxima); maxima = maxima + min(curTW) - 1;
%     [~,id] = max(x0(maxima));
%     tPeak = maxima(id);

    score = (x0(maxima)').*exp(-(maxima-mtPeak).^2/2/estPeakShift^2);
    [~,id] = max(score);
    tPeak = maxima(id); if(isempty(maxima)) tPeak = mtPeak; end;
    
    
    %% cut at min
    s0 = mean((x0(2:end)-x0(1:end-1)).^2)/2;
    stds(ii) = s0;
%     [vMin,tMin] = min(x0(curTW(1):tPeak));
%     tMin = curTW(1) + tMin - 1;
    tMin = tPeak;vMin = x0(tMin);
    ts = tMin;
    %% extend backward
    for t = tMin:-1:1
        if(x0(t)-vMin>=3*s0)
            break;
        else
           if(x0(t)<vMin) 
              vMin = x0(t) ;
              ts = t;
           end
        end
    end
%     x0(1:ts-1) = 0;
%     x0(tPeak+1:end) = 1;
%     vMax = x0(tPeak);
%     x0(ts:tPeak) = (x0(ts:tPeak)-vMin)/(vMax-vMin);
%     tEnd = max(tPeak,tEnd);
%     tst(ii,:) = x0;
%     ref(ii,:) = refBase;

    vMax = x0(tPeak);
    x0Org(1:ts-1) = vMin;
    x0Org(tPeak+1:end) = vMax;
    x0Org = x0Org - vMin;
    tEnd = max(tPeak,tEnd);
    tst(ii,:) = x0Org;
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
spMap1 = zeros(H0,W0);
% tEnd = min(T0,tEnd+5);
% tst = tst(:,1:tEnd);
% ref = ref(:,1:tEnd);
% refBase = refBase(:,1:tEnd);

for ii=1:numel(spLst)
    spMap1(spLst{ii}) = ii;
end
nSp = numel(spLst);
for ii=1:nSp
    sp0 = spLst{ii};
    [ih0,iw0] = ind2sub([H0,W0],sp0);
    neiPix = [];
    for jj=1:numel(dh)  % find neighbors in eight directions
        ih = max(min(H0,ih0+dh(jj)),1);
        iw = max(min(W0,iw0+dw(jj)),1);
        ihw = sub2ind([H0,W0],ih,iw);
        neiPix = union(neiPix,ihw);%;union(neiPix,setdiff(ihw,sp0));
    end
    neib0 = [];
    for k = (ii+1):nSp
        if(~isempty(intersect(neiPix,spLst{k})))
            neib0 = [neib0,k];
        end
    end
    
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
    [vMin,~] = min(x0Smo(t00:tPeak));
    thr = 0.3*(x0Smo(tPeak)-vMin);
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