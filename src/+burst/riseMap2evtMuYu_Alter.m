function [svLabel] = riseMap2evtMuYu_Alter(spLst,tDly,neibLst,sv_spLabels,maxRiseUnc,cDelay,sz)
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

% get local minimum as event starting point
seedLabel = [];
for ii=1:nSp
    neib0 = neibLst{ii};
    curDly = tDly(ii);
    neiDly = tDly(neib0);
    if(sum(curDly>neiDly)==0)
        seedLabel = [seedLabel;ii];
    end
end
sz = cellfun(@numel,spLst(seedLabel));
seedLabel = seedLabel(sz>30);

%% find outlier seeds
outlier = false(numel(seedLabel),1);
for i = 1:numel(seedLabel)
    neib0 = neibLst{seedLabel(i)};
    curDly = tDly(seedLabel(i));
    neiDly = tDly(neib0);
    dif = min(neiDly) - curDly;
    if(dif>cDelay)
        outlier(i) = true;
    end
end
trueSd = seedLabel(outlier);


%% remove outLier super voxel from neibLst
neibLst_rout = neibLst;
if(~isempty(trueSd))
    for i = 1:nSp
        neibLst_rout{i} = setdiff(neibLst_rout{i},trueSd);
    end
end

%% left part, remove weak local minimals
seedLabel = seedLabel(~outlier);
[~,seedOrd] = sort(riseX(seedLabel),'descend');
seedLabel = seedLabel(seedOrd);
LeftSeeds = seedLabel;
for i = 1:numel(seedLabel)
    curLabel = seedLabel(i);
    nonWeak = checkWhetherTrueSeed(curLabel,tDly,neibLst_rout,maxRiseUnc,LeftSeeds);
    if(~nonWeak)
        LeftSeeds = setdiff(LeftSeeds,curLabel);
    end
end
trueSd = [trueSd;LeftSeeds];

% dlyMap = inf(H,W);for i = 1:numel(spLst) dlyMap(spLst{i}) = tDly(i); end; figure;imagesc(dlyMap);colorbar;
% mmm = false(H,W); for i = 1:numel(trueSd) mmm(spLst{trueSd(i)}) = true; end;figure;imshow(mmm);
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
% seed inner link
mustLink = [];
rootlabels = zeros(numel(spLst),1);
for i = 1:numel(sv_spLabels)
    curSv_SpLabels = sv_spLabels{i};
    id1 = curSv_SpLabels(1);
    for k = 2:numel(curSv_SpLabels)
        id2 = curSv_SpLabels(k);
%         edgeMatrix = [edgeMatrix;[id1,id2,0]];
        id01 = UF_find(rootlabels,id1);
        id02 = UF_find(rootlabels,id2);
        rootlabels = linkNodes(id01,id02,rootlabels);
        mustLink = [mustLink;id1,id2];
    end
end

% % linking seeds for appplying MST
needRemoveEdges = [];
id1 = trueSd(1);
for k = 2:numel(trueSd)
    id2 = trueSd(k);
    id01 = UF_find(rootlabels,id1);
    id02 = UF_find(rootlabels,id2);
    rootlabels = linkNodes(id01,id02,rootlabels);
end
idLinkage = size(edgeMatrix,1)+[1:size(needRemoveEdges,1)];
edgeMatrix = [edgeMatrix;needRemoveEdges];

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
            newNeib = [newNeib;neibLst{neib0(i)}];
        end
        checked = [checked;neib0];
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


