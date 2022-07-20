function [spLst,cx,tDly,neiLst,sv_spLabels,rgt00,isFail] = spgtw_MuYuProject_Simply(...
    dF,seMap,seSel,smoBase,maxStp,cDelay,spSz,spT,superVoxels,opts)
% spgtw super pixel GTW 
% make one burst to super pixels and run gtw

if ~isfield(opts,'gtwGapSeedMin') || ~isfield(opts,'gtwGapSeedRatio')
    opts.gtwGapSeedRatio = 4;
    opts.gtwGapSeedMin = 5;
end

%% setting
[H,W,T] = size(dF);
isFail = 0;
maxStp = max(min(maxStp,round(T/2)),1);
% s00 = 1; % already normalized

%% super pixels
m0Msk = seMap==seSel;
validMap = sum(m0Msk,3)>0;
cc = bwconncomp(validMap); cc = cc.PixelIdxList;
sz = cellfun(@numel,cc); [~,id] = max(sz);
validMap = false(H,W); validMap(cc{id}) = true;

dFip = dF;
dFip(seMap~=seSel) = nan;
dFip = gtw.imputeMov_Fast2(dFip,validMap);

dFAvg = max(dFip,[],3);
dFAvg(isnan(dFAvg)) = 0;
% dFAvg(dFAvg<0) = 0;
% dFAvg = medfilt2(dFAvg);  

%% super pixel
spSzMinLimitation = spSz/4;
nSpEstimate = max(1,round(H*W/spSz));
[L] = superpixels(dFAvg,nSpEstimate);
L(~validMap) = 0;
spLst = label2idx(L);

%% get super pixel
if sum(validMap(:))< spSzMinLimitation || numel(spLst)<2
    spLst = {find(validMap>0)};
    cx = [];
    tDly = [];
    neiLst = [];
    rgt00 = 1:T;
    isFail = 1;
    sv_spLabels = [];
    return
end

%% connectivity
nSp = numel(spLst);
nAdd = nSp + 1;
for i = 1:nSp
    pix = spLst{i};
   if(~isempty(pix)) 
       cMap = false(H,W);
       cMap(pix) = true;
       cc = bwconncomp(cMap);
       if(cc.NumObjects>1)
           comps = cc.PixelIdxList;
           for j = 1:numel(comps)
               if(j==1)
                   spLst{i} = comps{j};
               else
                   spLst{nAdd} = comps{j};
                   L(comps{j}) = nAdd;
                   nAdd = nAdd+1;
               end
           end
       end
   end
end

%%
spLst = cell(0);
nSp = 0;
sv_spLabels = cell(numel(superVoxels),1);
for k = 1:numel(superVoxels)
    pix = superVoxels{k};
    [ih,iw,~] = ind2sub([H,W,T],pix);
    ihw = unique(sub2ind([H,W],ih,iw));
    L0 = zeros(H,W);
    L0(ihw) = L(ihw);
    curSp = mergeSmallSp(L0,spSzMinLimitation);
    spLabels = nSp + [1:numel(curSp)];
    nSp = nSp + numel(curSp);
    spLst(spLabels) = curSp;
    sv_spLabels{k} = spLabels;
end

%%
% signal part
idx0 = find(m0Msk>0);
[~,~,it0] = ind2sub(size(m0Msk),idx0);
rgt00 = max(min(it0)-5,1):min(max(it0)+5,T);
dat = imgaussfilt3(dFip(:,:,rgt00),[5,5,2]);
% dat(~m0Msk) = nan;
% dat = dat(:,:,rgt00);
datVec = reshape(dat,[],numel(rgt00));

fprintf('Node %d\n',numel(spLst));

%% get delays
% graph
dh = [-1 0 1 -1 0 1 -1 0 1];
dw = [-1 -1 -1 0 0 0 1 1 1];
vMap = seMap==seSel; vMap = vMap(:,:,rgt00); vMap = reshape(vMap,[],numel(rgt00));
tDly = zeros(nSp,1);
neiLst = cell(nSp,1);
thrVec = 0.8:0.02:0.95;

for i = 1:nSp
    pix = spLst{i};
    cx = mean(datVec(pix,:),1,'omitnan');
    cx = cx - min(cx);
    cx = cx/max(cx);
    tAch = zeros(numel(thrVec),1);
    for ii=1:numel(thrVec)
        tAch(ii) = find(cx>=thrVec(ii),1);
    end
%     cx(isnan(cx)) = 0; cx = imgaussfilt(cx,2);
%     validT = sum(vMap(pix,:)) == numel(pix);
    
%     [~,tPeak] = max(cx);
    tDly(i) = mean(tAch);
    
    [ih0,iw0] = ind2sub([H,W],pix);
    neiPix = [];
    for jj=1:numel(dh)  % find neighbors in eight directions
        ih = max(min(H,ih0+dh(jj)),1);
        iw = max(min(W,iw0+dw(jj)),1);
        ihw = sub2ind([H,W],ih,iw);
        neiPix = union(neiPix,ihw);%;union(neiPix,setdiff(ihw,sp0));
    end
    for k = (i+1):nSp
        if(~isempty(intersect(neiPix,spLst{k})))
            neiLst{i} = [neiLst{i};k];
            neiLst{k} = [neiLst{k};i];
        end
    end
end



end

function spLst = mergeSmallSp(L,spSzMinLimitation)
    spLst = label2idx(L);
    [H,W] = size(L);
    dh = [-1,-1,-1,0,0,1,1,1];
    dw = [-1,0,1,-1,1,-1,0,1];
    isolatedRegion = [];
    sz = cellfun(@numel,spLst);
    smallRegions = find(sz>0 & sz<spSzMinLimitation);
    % merge small pixels
    while(~isempty(smallRegions))
        for i = 1:numel(smallRegions)
            curLabel = smallRegions(i);
            pix = spLst{curLabel};
            curSz = numel(spLst{curLabel});
            if(curSz>0 && curSz<spSzMinLimitation)
                [ih,iw] = ind2sub([H,W],pix);
               for k = 1:numel(dw)
                   ih1 = min(H,max(1,ih+dh(k)));
                   iw1 = min(W,max(1,iw+dw(k)));
                   pixShift = sub2ind([H,W],ih1,iw1);
                   neiL = setdiff(L(pixShift),[0,curLabel]);
                   if(~isempty(neiL))
                      break; 
                   end
               end
               if(isempty(neiL))
                  isolatedRegion = [isolatedRegion,curLabel];
                  continue;
               end
               neiL = neiL(1);
               % update
               spLst{neiL} = [spLst{neiL};pix];
               L(pix) = neiL;
               spLst{curLabel} = [];
            end
        end
        sz = cellfun(@numel,spLst);
        smallRegions = find(sz>0 & sz<spSzMinLimitation);
        smallRegions = setdiff(smallRegions,isolatedRegion);
        if(isempty(smallRegions))
           break; 
        end
    end 
    sz = cellfun(@numel,spLst);
    spLst = spLst(sz>0);
end