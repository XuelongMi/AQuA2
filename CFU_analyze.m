
%% load res.mat
cfuInfo = res.cfuInfo1;
%%
pOut = 'D:\Test\'; %% tif folder
mkdir(pOut);

%% relation
tic;
nCFU = size(cfuInfo,1);
maxDist = 4;        % unfixed time window, pick the most significant one
relation = cell(nCFU*nCFU,1);
parfor k = 1:nCFU*nCFU
    i = floor((k-1)/nCFU)+1;
    j = k-(i-1)*nCFU;
    if(j<=i)
        relation{k} = [];
    else
        seq1 = cfuInfo{i,4};
        seq2 = cfuInfo{j,4};
        [pvalue1,ds1,distribution1] = cfu.calDependency(seq1, seq2,0:maxDist); % condition is the first variable, occurrence is the second.
        [pvalue2,ds2,distribution2] = cfu.calDependency(seq2, seq1,0:maxDist); % condition is the first variable, occurrence is the second.
        delay = nan;
        if(pvalue1<pvalue2)
            pvalue = pvalue1;
            if(~isempty(distribution1))
                delay = distribution1(:,1)'*distribution1(:,2)/sum(distribution2(:,2));
            end
        else
            pvalue = pvalue2;
            if(~isempty(distribution2))
                delay = -distribution2(:,1)'*distribution2(:,2)/sum(distribution2(:,2));
            end
        end
        relation{k} = [i,j,pvalue,delay];
    end
end

valid = false(nCFU*nCFU,1);
parfor k = 1:nCFU*nCFU
    if(~isempty(relation{k})&& relation{k}(3)<1e-3)
        valid(k) = true;
    end
end

relation = relation(valid);

save([pOut,'cfuInfo.mat'],'cfuInfo','-v7.3');
save([pOut,'cfuInfoRelation.mat'],'relation','-v7.3');
toc;

%% clustering
thr = 1e-2;
cfuNumThr = 3;

% filter
nPair = 0;
cnt = 0;
relation = cell2mat(relation);
select = relation(:,3)<thr;
relation = relation(select,:);

% preprocess
select = relation(:,4)<0;
relation(select,[1,2]) = relation(select,[2,1]);
relation(select,4) = -relation(select,4);

%
% sort
[~,id] = sort(relation(:,3));
relation = relation(id,:);
groupLabels = 1:nCFU;
addedEdge = nan(nCFU,1);
% clustering
for k = 1:size(relation,1)
    id01 = relation(k,1);
    id02 = relation(k,2);
    id1 = cfu.findRootLabel(groupLabels,id01);
    id2 = cfu.findRootLabel(groupLabels,id02);
    groupLabels(id1) = min(id1,id2);
    groupLabels(id2) = min(id1,id2);
end

for k = 1:nCFU
    root = cfu.findRootLabel(groupLabels,k);
    groupLabels(k) = root;
end

% sort
cc = label2idx(groupLabels);
id = cellfun(@numel,cc)>cfuNumThr;
cc = cc(id);
[~,id] = sort(cellfun(@numel,cc),'descend');
cc = cc(id);

% groupInfo
visited = false(nCFU,1);
groupInfo = cell(numel(cc),4);
delays = nan(nCFU,1);
addedPvalue = nan(nCFU,1);
for k = 1:numel(cc)
   groupInfo{k,1} = k;
   labels = cc{k};
   p_values = zeros(numel(labels),1);
   for i = 1:numel(labels)
      curLabel = labels(i) ;
      select = relation(:,1)==curLabel | relation(:,2)==curLabel;
      p_values(i) = mean(relation(select,3));
   end
   [~,id] = min(p_values);
   curLabel = labels(id);
   delays(curLabel) = 0;
   alreadyIncluded = curLabel;
   addLst = curLabel;
   while(~isempty(addLst))
       curLabel = addLst(1);
       %
       id =  find(relation(:,1) == curLabel);
       for i = 1:numel(id)
           id2 = relation(id(i),2);
           if(~visited(id2))
               delays(id2) = relation(id(i),4) + delays(curLabel);
               addedPvalue(id2) = relation(id(i),3);
               addLst = [addLst,id2];
               visited(id2) = true;
           end
       end
       
       id =  find(relation(:,2) == curLabel);
       for i = 1:numel(id)
           id1 = relation(id(i),1);
           if(~visited(id1))
               delays(id1) = - relation(id(i),4) + delays(curLabel);
               addedPvalue(id1) = relation(id(i),3);
               addLst = [addLst,id1];
               visited(id1) = true;
           end
       end
       addLst = addLst(2:end);
   end
   
   groupInfo{k,2} = labels;
   groupInfo{k,3} = delays(labels);
   groupInfo{k,4} = addedPvalue(labels);
end
save([pOut,'cfuGroupInfo.mat'],'groupInfo','-v7.3');

%% show groups with delay info
close all;
datPro = mean(res.datOrg1,3);
datPro = datPro - mean(datPro(:));
datPro = datPro / max(datPro(:));
[H,W,nSlice] = size(datPro);
colorMap = parula(256);
nGroup = size(groupInfo,1);
minDelayToShows = zeros(nGroup,1);
maxDelayToShows = zeros(nGroup,1);
for kk = 1:nGroup
    groupMember = groupInfo{kk,2};
    relativeDelay = groupInfo{kk,3};
    
    minDelayToShows(kk) = min(relativeDelay);
    maxDelayToShows(kk) = max(relativeDelay);
    
    map = zeros(H,W,3,nSlice);
    for k = 1:nSlice
        map(:,:,:,k) = cat(3,datPro(:,:,k),datPro(:,:,k),datPro(:,:,k));
    end
    
    colorId = round((relativeDelay-minDelayToShows(kk))/(max(relativeDelay)-minDelayToShows(kk))*255+1);
    colorId(isnan(colorId)) = 128;
    for k = 1:numel(groupMember)
        curID = groupMember(k);
        curZ = cfuInfo{curID,2};
        rcolor = colorMap(colorId(k),1);
        gcolor = colorMap(colorId(k),2);
        bcolor = colorMap(colorId(k),3);
        curCFU = cfuInfo{curID,3};
        curCFU(isnan(curCFU))=0;
        map(:,:,1,curZ) = map(:,:,1,curZ) + rcolor*(curCFU);
        map(:,:,2,curZ) = map(:,:,2,curZ) + gcolor*(curCFU);
        map(:,:,3,curZ) = map(:,:,3,curZ) + bcolor*(curCFU);
    end
    zzshow(map);
end
