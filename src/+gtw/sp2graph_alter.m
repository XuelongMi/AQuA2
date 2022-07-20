function [ref,tst,refBase,s,t,idxGood] = sp2graph_alter(df0ip,vMap0,spLst,seedIn,gapSeedHW)
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
df00Vec = reshape(df0ip(rgh,rgw,:),[],T0);
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
% df0VecSmo = reshape(df0ipSmo,[],T0);

nSp = numel(spLst);
tst = zeros(nSp,T0);
ref = zeros(nSp,T0);
for ii=1:numel(spLst)
    sp0 = spLst{ii};
    %     if numel(sp0)>2
    % scale and baseline
%     tst0smo = nanmean(df0VecSmo(sp0,:),1);
%     tst0smo = tst0smo - min(tst0smo);
%     k0 = max(tst0smo)/max(refBase);
    
    tst0 = nanmean(df0Vec(sp0,:),1);
    tst0g = imgaussfilt(tst0,1);
%     tst0 = tst0g - min(tst0g);
    [tst0,tPeak] = balanceCurve(tst0g);
    tEnd = max(tPeak,tEnd);
    ref0 = refBase*max(tst0)/max(refBase);
    tst(ii,:) = tst0;
    ref(ii,:) = ref0;
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
tEnd = min(T0,tEnd+5);
tst = tst(:,1:tEnd);
ref = ref(:,1:tEnd);
refBase = refBase(:,1:tEnd);

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
% function curve = balanceCurve(curve)
%     s0 = mean((curve(2:end)-curve(1:end-1)).^2)/2;
%     v1 = curve(1);
%     v2 = curve(end);
%     t1 = find(curve>v1+3*s0,1);
%     t2 = find(curve>v2+3*s0,1,'last');
%     curve(1:t1-1) = 0;
%     curve(t2+1:end) = 0;
%     if(t2~=t1)
%         curve(t1:t2) = curve(t1:t2) - (curve(t2)-curve(t1))/(t2-t1)*([t1:t2]-t1) - curve(t1);
%     end
% end
function [curve,tPeak] = balanceCurve(curve)
    [vP,tPeak] = max(curve);
    [vM,tMin] = min(curve(1:tPeak));
    curve(1:tMin) = vM;
    curve(tPeak:end) = vP;
    curve = curve - vM;
end