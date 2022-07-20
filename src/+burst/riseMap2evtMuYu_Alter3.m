function [svLabel] = riseMap2evtMuYu_Alter3(spLst,tDly,neibLst,sv_spLabels,maxRiseUnc,cDelay,sz,spSz,riseDur)
% cluster super voxels to events
% spLst: super voxel list
% tDly:  rising time of super pixels
% neibLst: the neighbor relationship between super pixels
% sv_spLabels: the relationship between super pixels and super voxels
% maxRiseUnc:  the slack, or the uncertainty for determining delay time
% cDelay: the rising time difference for determining outliers



%%
H = sz(1);W = sz(2);
nSp = numel(spLst);
riseX = tDly;

%% sv,sp relation
sp_svLabel = zeros(numel(nSp),1);
for i = 1:numel(sv_spLabels)
    curLabels = sv_spLabels{i};
    sp_svLabel(curLabels) = i;
end

% %% find outlier seeds
OutSeeds = [];
for i = 1:nSp
    neib0 = neibLst{i};
    curDly = tDly(i);
    neiDly = tDly(neib0);
    dif = min(neiDly) - curDly;
    if(isempty(neib0))
        OutSeeds = [OutSeeds;i];
    end
end

%% merge seeds if belong to same super voxels
OutSeeds = mergeSeeds(OutSeeds,sp_svLabel);

%% remove outLier super voxel from neibLst, isolate outlier
removeNei = [];
for i = 1:numel(OutSeeds)
    cSv = sp_svLabel(OutSeeds(i));
    curLabel = sv_spLabels{cSv};
    removeNei = [removeNei,curLabel];
end
removeLabels = false(nSp,1);
removeLabels(removeNei) = true;
if(~isempty(OutSeeds))
    for i = 1:nSp
        if(~removeLabels(i))
            neibLst{i} = setdiff(neibLst{i},removeNei);
        else
            cSv = sp_svLabel(i);
            neibLst{i} = sv_spLabels{cSv};
        end
    end
end

% get local minimum as event starting point, except outliers
seedLabel = [];
for ii=1:nSp
    if(~removeLabels(ii))
        neib0 = neibLst{ii};
        curDly = tDly(ii);
        neiDly = tDly(neib0);
        dif = min(neiDly) - curDly;
        if(sum(curDly>neiDly)==0 && (dif<=cDelay || numel(spLst{ii})>4*spSz))
            seedLabel = [seedLabel;ii];
        end
    end
end
sz = cellfun(@numel,spLst(seedLabel));
seedLabel = seedLabel(sz>30);

%% left part, remove weak local minimals
[~,seedOrd] = sort(riseX(seedLabel),'descend');
seedLabel = seedLabel(seedOrd);
LeftSeeds = seedLabel;
for i = 1:numel(seedLabel)
    curLabel = seedLabel(i);
%     maxRiseUnc0 = max(riseDur(i)*0.1,maxRiseUnc);
    maxRiseUnc0 = maxRiseUnc;
    nonWeak = checkWhetherTrueSeed(curLabel,tDly,neibLst,maxRiseUnc0,LeftSeeds);
    if(~nonWeak)
        LeftSeeds = setdiff(LeftSeeds,curLabel);
    end
end
LeftSeeds = mergeSeeds(LeftSeeds,sp_svLabel);
trueSd = [OutSeeds;LeftSeeds];

% dlyMap = nan(H,W);for i = 1:numel(spLst) dlyMap(spLst{i}) = tDly(i); end; figure;imagesc(dlyMap);colorbar;
% mmm = zeros(H,W); for i = 1:numel(trueSd) mmm(spLst{trueSd(i)}) = trueSd(i); end;zzshow(mmm);

%% make sure each component has one seed
ccGraph = zeros(nSp,1);
for i = 1:nSp
    neib0 = neibLst{i};
    id1 = UF_find(ccGraph,i);
    for j = 1:numel(neib0)
        id2 = neib0(j);
        id2 = UF_find(ccGraph,id2);
        if(id1~=id2)
            ccGraph = linkNodes(id1,id2,ccGraph);
        end
    end
end
for i = 1:nSp
    ccGraph(i) = UF_find(ccGraph,i);
end
seLb = false(nSp,1);    seLb(trueSd) = true;
[uccGraph] = unique(ccGraph);
for i = 1:numel(uccGraph)
    curCom = find(ccGraph==uccGraph(i));
    if(sum(seLb(curCom))==0)
        [~,id] = min(riseX(curCom));
        trueSd = [OutSeeds;curCom(id)];
    end
end
%% construct graph
% super pixel relationship
edgeMatrix = [];
for i = 1:nSp
    neib0 = neibLst{i};
    neib0 = neib0(neib0>i);
    if(size(neib0,1)==1) neib0 = neib0';end
    dif = abs(riseX(neib0)-riseX(i)) + 1e-6;
    addMatrix = [ones(numel(neib0),1)*i,neib0,dif];
    edgeMatrix = [edgeMatrix;addMatrix];
end

% linkage for each super voxel
mustLink = [];
rootlabels = zeros(numel(spLst),1);
for i = 1:numel(sv_spLabels)
    curSv_SpLabels = sv_spLabels{i};
    id1 = curSv_SpLabels(1);
    for k = 2:numel(curSv_SpLabels)
        id2 = curSv_SpLabels(k);
        id01 = UF_find(rootlabels,id1);
        id02 = UF_find(rootlabels,id2);
        rootlabels = linkNodes(id01,id02,rootlabels);
        mustLink = [mustLink;id1,id2];
    end
end

% % linking seeds for appplying MST
id1 = trueSd(1);
for k = 2:numel(trueSd)
    id2 = trueSd(k);
    id01 = UF_find(rootlabels,id1);
    id02 = UF_find(rootlabels,id2);
    rootlabels = linkNodes(id01,id02,rootlabels);
end

%% MST
validEdges = false(size(edgeMatrix,1),1);
if(~isempty(edgeMatrix))
    [~,id] = sort(edgeMatrix(:,3));
    edgeMatrix = edgeMatrix(id,:);
end

for i = 1:size(edgeMatrix,1)
    id1 = edgeMatrix(i,1);
    id2 = edgeMatrix(i,2);
    id1 = UF_find(rootlabels,id1);
    id2 = UF_find(rootlabels,id2);
    if(id1~=id2)
        rootlabels = linkNodes(id1,id2,rootlabels);
        validEdges(i) = true;
    end
end
%% remove seed linkage then link again
if(~isempty(edgeMatrix))
    mustLink = [mustLink;edgeMatrix(validEdges,1:2)];
end
rootlabels2 = zeros(nSp,1);
for i = 1:size(mustLink,1)
    id1 = mustLink(i,1);
    id2 = mustLink(i,2);
    id1 = UF_find(rootlabels2,id1);
    id2 = UF_find(rootlabels2,id2);
    rootlabels2 = linkNodes(id1,id2,rootlabels2);
end
for i = 1:nSp
    rootlabels2(i) = UF_find(rootlabels2,i);
end

%% update
[~,ia,ic] = unique(rootlabels2);
sp_evt = zeros(nSp,1);
for i = 1:nSp
   seID = ic(i);
   sp_evt(i) = seID;
end
svLabel = zeros(numel(sv_spLabels),1);
for i = 1:numel(sv_spLabels)
    labels = sv_spLabels{i};
    svLabel(i) = sp_evt(labels(1));
end
    
end
function root = UF_find(labels,id)
   
   if(labels(id)==0 || labels(id)==id)
        root = id;
   else
        root = UF_find(labels,labels(id));
   end
    
end
function nonWeak = checkWhetherTrueSeed(curLabel,tDly,neibLst,maxRiseUnc,LeftSeeds)
    nonWeak = false;
    % BFS search
    tDlyThr = tDly(curLabel) + maxRiseUnc;
    neib0 = neibLst{curLabel};
    neib0 = neib0(tDly(neib0)<=tDlyThr);
    checked = [curLabel];
    while(~isempty(neib0))
        if(numel(intersect(LeftSeeds,neib0))>0)
           return; 
        end
        newNeib = [];
        for i = 1:numel(neib0)
            newNeib = [newNeib,neibLst{neib0(i)}];
        end
        checked = [checked,neib0];
        newNeib = setdiff(newNeib,checked);
        newNeib = newNeib(tDly(newNeib)<=tDlyThr);
        neib0 = newNeib;
    end
    nonWeak = true;
end
function labels = linkNodes(id1,id2,labels)
    root = min(id1,id2);
    labels(id1) = root;
    labels(id2) = root;
end
function seedLst = mergeSeeds(seedLst,sp_svLabels)
    realSeeds = false(numel(sp_svLabels),1);
    nSv = max(sp_svLabels);
    detected = false(nSv,1);
    for i = 1:numel(seedLst)
        cSv = sp_svLabels(seedLst(i));
        if(~detected(cSv))
            detected(cSv) = true;
            realSeeds(i) = true;
        end
    end
    seedLst = seedLst(realSeeds);
end

