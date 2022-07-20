function [svLabel,delayMap] = delayMap2evt(spLst,tDly,majorityMap,major0,spSz,opts)
% cluster super voxels to events
% spLst: super voxel list
% tDly:  rising time of super pixels
% neibLst: the neighbor relationship between super pixels
% sv_spLabels: the relationship between super pixels and super voxels
% maxRiseUnc:  the slack, or the uncertainty for determining delay time
% cDelay: the rising time difference for determining outliers

[H,W] = size(majorityMap);
delayMap = nan(H,W);
for i = 1:numel(spLst)
   delayMap(spLst{i}) = tDly(i);
end

dh = [-1,-1,-1,0,0,1,1,1];
dw = [-1,0,1,-1,1,-1,0,1];
minDly = min(tDly);
maxDly = max(tDly);
thrs = minDly:(maxDly-minDly)/10:maxDly;
seedMap = zeros(H,W);
nSeed = 0;
for k = 1:numel(thrs)    
    thr = thrs(k);
    candidateRegions = bwconncomp(delayMap<thr);
    candidateRegions = candidateRegions.PixelIdxList;
    sz = cellfun(@numel,candidateRegions);
    candidateRegions = candidateRegions(sz>3*spSz);
    for i = 1:numel(candidateRegions)
        pix = candidateRegions{i};
        seedsInRegion = setdiff(seedMap(pix),0);
        if(numel(seedsInRegion)==1)
            seedMap(pix) = seedsInRegion(1);
        elseif(numel(seedsInRegion)==0)
            newAdd = pix;
            neighbor = [];
            pixGrow = pix;
            round = 0;
            while(round<100 && numel(newAdd)>0 && numel(neighbor)<numel(pix))
                [ih0,iw0] = ind2sub([H,W],newAdd);
                newAdd = [];
                for ii = 1:numel(dh)
                    ih = min(max(1,ih0 + dh(ii)),H);
                    iw = min(max(1,iw0 + dw(ii)),W);
                    newAdd = [newAdd;sub2ind([H,W],ih,iw)];
                end
                newAdd = setdiff(newAdd,pixGrow);
                newAdd = newAdd(~isnan(delayMap(newAdd)));
                neighbor = [neighbor;newAdd];
                pixGrow = [pixGrow;newAdd];
                round = round + 1;
            end
            if(mean(delayMap(neighbor)) - nanmean(delayMap(pix))>opts.cDelay)
                nSeed = nSeed + 1;
                seedMap(pix) = nSeed;
            end
        end
    end
end

majorLst = label2idx(majorityMap);
svLabel = zeros(numel(majorLst),1);
if(max(seedMap(:))==0)
    svLabel(:) = 1;
    return;
end

scoreMap0 = delayMap;
sdLst = label2idx(seedMap);

% deal with background
SE = strel('disk',2);
BW = true(size(scoreMap0));
BW(majorityMap>-1) = false;
BW2 = imerode(BW,SE);
scoreMap0(BW) =  max(tDly);
scoreMap0(BW2) = 0;
scoreMap1 = imimposemin(scoreMap0,seedMap>0|BW2,8);
MapOut = watershed(scoreMap1,8);
% update
MapOut(BW) = 0;
waterLst = label2idx(MapOut);
for ii = 1:numel(sdLst)
   target = sdLst{ii}(1);
   waterLabel = MapOut(target);
   seedMap(waterLst{waterLabel}) = ii;
end

for i = 1:numel(major0)
   ihw = major0{i}.ihw;
   svLabel(i) = mode(seedMap(ihw));
end

lst = label2idx(svLabel);
sz = cellfun(@numel,lst);
lst = lst(sz>0);
for i = 1:numel(lst)
    svLabel(lst{i}) = i;
end

end

