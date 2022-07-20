function [cfu_pre] = CFU_tmp_function(evts,sz,ff)
    H = sz(1);
    W = sz(2);
    T = sz(3);
    evtMap = zeros([H,W,T],'uint16');
    for i = 1:numel(evts)
       evtMap(evts{i}) = i;
    end
    evtMap = reshape(evtMap,[],T);

    %% tracking setting
    nNode = numel(evts);
    evtIhw = cell(nNode,1);
    weightedIhw = cell(nNode,1);
    maxCounts = zeros(nNode,1);
    parfor i = 1:numel(evtIhw)
       [ih,iw,it]  = ind2sub([H,W,T],evts{i});
       [count,ihw] = groupcounts(sub2ind([H,W],ih,iw));
       evtIhw{i} = ihw;
       weightedIhw{i} = count/max(count);
       maxCounts(i) = max(count);
    end
    nNode = numel(evts);
    %% linkage
    
    s_t0 = zeros(nNode*nNode,3);
    overlapLst = cell(nNode,1); % to avoid circle loop
    nPair = 0;
    tic;
    for i = 1:nNode
        waitbar(0.2 + i/nNode*0.8,ff);
        label1 = i;
        ihw1 = evtIhw{i};
        weightedIhw1 = weightedIhw{i};
        
        possibleCandidates = setdiff(evtMap(ihw1,:),[0,i]);
        overlapLst{i} = possibleCandidates;
        possibleCandidates = possibleCandidates(possibleCandidates>i);
        nPC = numel(possibleCandidates);
        if(size(possibleCandidates,1)==1)
            possibleCandidates = possibleCandidates';
        end
        overlap = zeros(nPC,1);
        %% accelerate
        for j = 1:nPC
            ihw2 = evtIhw{possibleCandidates(j)};
            weightedIhw2 = weightedIhw{possibleCandidates(j)};
            Un = union(ihw1,ihw2);
            w1 = zeros(numel(Un),1);
            w1(ismember(Un,ihw1)) = weightedIhw1;
            w2 = zeros(numel(Un),1);
            w2(ismember(Un,ihw2)) = weightedIhw2;
            overlap(j) = sum(min(w1,w2))./sum(max(w1,w2));
        end
        range = nPair + [1:nPC];
        s_t0(range,:) = [ones(nPC,1)*label1,double(possibleCandidates),overlap];
        nPair = nPair + nPC;
    end
    s_t0 = s_t0(1:nPair,:);
    toc;    
    cfu_pre.s_t0 = s_t0;
    cfu_pre.weightedIhw = weightedIhw;
    cfu_pre.maxCounts = maxCounts;
    cfu_pre.overlapLst = overlapLst;
    cfu_pre.evtIhw = evtIhw;
end