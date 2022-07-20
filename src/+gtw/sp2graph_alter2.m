function [ref,tst,refBase,s,t,idxGood,orgCurves] = sp2graph_alter2(rDF,df0ip,m0Msk,vMap0,spLst,seedIn,gapSeedHW)
% sp2graph convert super pixels to curves and graph nodes for GTW
% Input super pixels are not perfect:
% Empty or too small
% Bad corresponding curves

[H0,W0,T0] = size(df0ip);

% test and refererence curves
% optionally, ignore signals after arriving peak
% df0ipSmo = imgaussfilt3(df0ip,[1 1 1]);

[ih,iw] = ind2sub([H0,W0],seedIn);
rgh = max(ih-gapSeedHW,1):min(ih+gapSeedHW,H0);
rgw = max(iw-gapSeedHW,1):min(iw+gapSeedHW,W0);
df00Vec = reshape(rDF(rgh,rgw,:),[],T0);
vm00 = vMap0(rgh,rgw);
df00Vec = df00Vec(vm00>0,:);
refBase = nanmean(df00Vec,1);
refBase = imgaussfilt(refBase,1);
% refBase = refBase - nanmin(refBase);


r1 = refBase;
[~,ix] = max(r1);
bw = r1*0;
bw(ix) = 1;
r2 = -imimposemin(-r1,bw);
r2(ix) = r1(ix);
r2(isinf(r2)) = nan;
refBase = r2;
[refBase,tEnd] = balanceCurve(refBase);

df0Vec = reshape(df0ip,[],T0);
m0Msk = reshape(m0Msk,[],T0);
% df0VecSmo = reshape(df0ipSmo,[],T0);

nSp = numel(spLst);
tst = zeros(nSp,T0);
ref = zeros(nSp,T0);
orgCurves = zeros(nSp,T0);
for ii=1:numel(spLst)
    sp0 = spLst{ii};
    valid = find(sum(m0Msk(sp0,:))>0);
    t0 = valid(1);
    x0 = mean(df0Vec(sp0,:),1);
    x0 = imgaussfilt(x0,2);
    [vMax,tPeak] = max(x0(valid));
    orgCurves(ii,:) = x0;
%     tPeak = valid(tPeak);
    
    curTW = max(min(valid)-2,1):min(T0,max(valid)+2);
    left = false(numel(curTW),1);  left(2:end) = x0(min(curTW)+1:max(curTW))>=x0(min(curTW):max(curTW)-1);
    right = false(numel(curTW),1); right(1:end-1) = x0(min(curTW):max(curTW)-1)>=x0(min(curTW)+1:max(curTW));
    maxima = left & right; maxima = find(maxima); maxima = maxima + min(curTW) - 1;
    [~,id] = max(x0(maxima));
    tPeak = maxima(id);
    
    %% cut at min
    s0 = mean((x0(2:end)-x0(1:end-1)).^2)/2;
    [vMin,tMin] = min(x0(t0:tPeak));
    tMin = t0 + tMin - 1;
%     tMin = tPeak;vMin = x0(tMin);
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
    x0(1:ts-1) = 0;
    x0(tPeak+1:end) = 1;
    x0(ts:tPeak) = (x0(ts:tPeak)-vMin)/(vMax-vMin);
    tEnd = max(tPeak,tEnd);
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