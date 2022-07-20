function [evtMap] = riseMap2evtMuYu(spLst,dlyMap,distMat,maxRiseUnc,cDelay)

nSp = numel(spLst);
[H,W] = size(dlyMap);

if(nSp==0)
    evtMemC = [];
    evtMemCMap = [];
   return; 
end

spMap = zeros(H,W);
riseX = nan(nSp,1);
for nn=1:nSp
    spMap(spLst{nn}) = nn;
    riseX(nn) = dlyMap(spLst{nn}(1));
end

% get local minimum as event starting point
% small area starting point not trustable, merged by its neighbors
initAreaThr = 30;
dlyMap1 = dlyMap;
mskDly = ~isinf(dlyMap1);
dw = [-1,-1,-1,0,0,1,1,1];
dh = [-1,0,1,-1,1,-1,0,1];
for ii=1:1000
    lm = imregionalmin(round(dlyMap1*100));
    lm = lm.*mskDly;
    if max(lm(:))==0
        lm = mskDly;
    end
    cc = bwconncomp(lm);
    allOK = 1;
    for jj=1:cc.NumObjects
        pix0 = cc.PixelIdxList{jj};
        if numel(pix0)<initAreaThr
            sp00 = spMap(pix0);
            sp00 = unique(sp00);
            sp00 = sp00(sp00>0);
            xNeib = [];
            [ih,iw] = ind2sub([H,W],pix0);
            for kk = 1:numel(dw)
                ih1 = max(1,min(ih+dh(kk),H));
                iw1 = max(1,min(iw+dw(kk),W));
                pix_change = sub2ind([H,W],ih1,iw1);
                xNeib = union(xNeib,setdiff(spMap(pix_change),[0,sp00]));
            end
            if ~isempty(xNeib)
                allOK = 0;
                tNew = min(riseX(xNeib));
                riseX(sp00) = tNew;
                dlyMap1(pix0) = tNew;
            end
        end
    end
    if allOK==1
        break
    end
end

sdLst = cc.PixelIdxList;
%% find outlier seeds
outlier = false(numel(sdLst),1);
trueSd = cell(0);
nSd = 1;
for i = 1:numel(sdLst)
    pix0 = sdLst{i};
    [ih,iw] = ind2sub([H,W],pix0);
    pixN = [];
    for kk = 1:numel(dw)
        ih1 = max(1,min(ih+dh(kk),H));
        iw1 = max(1,min(iw+dw(kk),W));
        pix_change = sub2ind([H,W],ih1,iw1);
        pixN = union(pixN,setdiff(pix_change,pix0));
    end
    dif = min(dlyMap1(pixN)) - dlyMap1(pix0(1));
    if(dif>cDelay)
        trueSd{nSd} = pix0;
        nSd = nSd + 1;
        outlier(i) = true;
    end
end

%% left part, remove weak local minimals
sdLst = sdLst(~outlier);
if(~isempty(sdLst))
    riseSd = zeros(numel(sdLst),1);
    for i = 1:numel(sdLst)
        riseSd(i) = dlyMap1(sdLst{i}(1));
    end
    [riseSd,seedOrd] = sort(riseSd,'descend');
    sdLst = sdLst(seedOrd);

    sdMap = zeros(H,W);
    for i = 1:numel(sdLst)
       sdMap(sdLst{i}) = i;
    end

    for i = 1:numel(sdLst)
        pix0 = sdLst{i};
        cc = bwconncomp(dlyMap1<=riseSd(i)+maxRiseUnc);
        cc = cc.PixelIdxList;
        curRegion = [];
        for j = 1:numel(cc)
           if(ismember(pix0(1),cc{j})) 
               curRegion = cc{j};
               break;
           end
        end
        sdLabels = setdiff(sdMap(curRegion),0);
        if(numel(sdLabels)==1)
            trueSd{nSd} = pix0;
            nSd = nSd + 1;
        else
            sdMap(pix0) = 0;
        end
    end
else    % except outlier, no seeds
    dlyMap00 = dlyMap1;
    for i = 1:numel(trueSd)
        dlyMap00(trueSd{i}) = inf;
    end
    pix = find(~isinf(dlyMap00));
    trueSd{nSd} = pix;
    nSd = nSd + 1;
end

%%
% sdMap = zeros(H,W);
% for i = 1:numel(trueSd)
%     sdMap(trueSd{i}) = i;
% end
% spEvt = zeros(numel(spLst),1);
% for i = 1:numel(spLst)
%    pix0 = spLst{i} ;
%    curLabel = mode(sdMap(pix0));
%    spEvt(i) = curLabel;
% end
% 
% distMat = nan(nSp,nSp);
% for i = 1:numel(spLst)
%     pix0 = spLst{i};
%     [ih,iw] = ind2sub([H,W],pix0);
%     xNeib = [];
%     for kk = 1:numel(dw)
%         ih1 = max(1,min(ih+dh(kk),H));
%         iw1 = max(1,min(iw+dw(kk),W));
%         pix_change = sub2ind([H,W],ih1,iw1);
%         xNeib = union(xNeib,setdiff(spMap(pix_change),[0]));
%     end
%     distMat(i,xNeib) = abs(riseX(i)-riseX(xNeib));
%     distMat(xNeib,i) = abs(riseX(i)-riseX(xNeib));
% end
% 
% 
% spEvt = burst.evtGrowLm(spEvt,distMat,riseX,spMap);
% 
% evtMap = zeros(H,W);
% % gather events
% pixLst = label2idx(spMap);
% evt0 = unique(spEvt);
% evt0 = evt0(evt0>0);
% nEvt0 = numel(evt0);
% for ii=1:nEvt0
%     pix0 = pixLst(spEvt==evt0(ii));
%     for jj=1:numel(pix0)
%         pix00 = pix0{jj};
%         evtMap(pix00) = ii;
%     end
% end
%  figure;imagesc(evtMap);colorbar


% %% construct graph, seed relationship
sdMap = zeros(H,W);
for i = 1:numel(trueSd)
    sdMap(trueSd{i}) = i;
end
sdLabels = cell(numel(trueSd),1);
for i = 1:numel(spLst)
   pix0 = spLst{i} ;
   curLabel = mode(sdMap(pix0));
   if(curLabel>0)
       sdLabels{curLabel} = [sdLabels{curLabel},i];
   end
end

% sp relation
edgeMatrix = [];
for i = 1:numel(spLst)
    pix0 = spLst{i};
    [ih,iw] = ind2sub([H,W],pix0);
    xNeib = [];
    for kk = 1:numel(dw)
        ih1 = max(1,min(ih+dh(kk),H));
        iw1 = max(1,min(iw+dw(kk),W));
        pix_change = sub2ind([H,W],ih1,iw1);
        xNeib = union(xNeib,setdiff(spMap(pix_change),[0,i]));
    end
    xNeib = xNeib(xNeib>i);
    if(size(xNeib,1)==1) xNeib = xNeib';end
    dif = abs(riseX(xNeib)-riseX(i)) + 1e-6;
    addMatrix = [ones(numel(xNeib),1)*i,xNeib,dif];
    edgeMatrix = [edgeMatrix;addMatrix];
end

sz = cellfun(@numel,sdLabels);
sdLabels = sdLabels(sz>0);
% seed inner link
for i = 1:numel(sdLabels)
    curSeedLabels = sdLabels{i};
    id1 = curSeedLabels(1);
    for k = 2:numel(curSeedLabels)
        id2 = curSeedLabels(k);
        edgeMatrix = [edgeMatrix;[id1,id2,0]];
    end
end

% seed intra link
needRemoveEdges = [];
id1 = sdLabels{1}(1);
for k = 2:numel(sdLabels)
    id2 = sdLabels{k}(1);
    needRemoveEdges = [needRemoveEdges;[id1,id2,0]];
end
idLinkage = size(edgeMatrix,1)+[1:size(needRemoveEdges,1)];
edgeMatrix = [edgeMatrix;needRemoveEdges];

%% MST
validEdges = false(size(edgeMatrix,1),1);
[~,id] = sort(edgeMatrix(:,3));
edgeMatrix = edgeMatrix(id,:);
rootlabels = zeros(numel(spLst),1);
for i = 1:size(edgeMatrix,1)
    id1 = edgeMatrix(i,1);
    id2 = edgeMatrix(i,2);
    id1 = UF_find(rootlabels,id1);
    id2 = UF_find(rootlabels,id2);
    
    if(id1~=id2)
        root = min(id1,id2);
        rootlabels(id1) = root;
        rootlabels(id2) = root;
        validEdges(i) = true;
    end
end
%% remove seed linkage
validEdges(ismember(id,idLinkage)) = false;
rootlabels2 = zeros(nSp,1);
for i = 1:size(edgeMatrix,1)
   if(validEdges(i)) 
       id1 = edgeMatrix(i,1);
       id2 = edgeMatrix(i,2);
       id1 = UF_find(rootlabels2,id1);
       id2 = UF_find(rootlabels2,id2);
       root = min(id1,id2);
       rootlabels2(id1) = root;
       rootlabels2(id2) = root;
   end
end
for i = 1:nSp
    rootlabels2(i) = UF_find(rootlabels2,i);
end

%% update
evtMap = zeros(H,W);
[~,ia,ic] = unique(rootlabels2);
for i = 1:nSp
   seID = ic(i);
   evtMap(spLst{i}) = seID;
end    
    
end
function root = UF_find(labels,id)
   
   if(labels(id)==0 || labels(id)==id)
        root = id;
   else
        root = UF_find(labels,labels(id));
   end
    
end



