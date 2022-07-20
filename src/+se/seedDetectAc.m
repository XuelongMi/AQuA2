function [sdLst,Map] = seedDetectAc(dF,dFResize,activeMap,opts,ff)
%% Detect seeds
    Thrs = opts.maxdF1:-opts.ratio:opts.thrARScl;
    [H,W,T] = size(dF);
    scaleRatios = opts.scaleRatios;
    Map = zeros(size(dF),'uint16');
    nEvt = 1;
    
%     time = 0;
    for k = 1:numel(Thrs)
        if exist('ff','var')
            waitbar(0.1 + 0.3*k/numel(Thrs),ff);
        end
        
        %% Find candidates
        curThr = Thrs(k);
        selectMap = (dF>curThr) & activeMap;
        curRegions = bwconncomp(selectMap);
        curRegions = curRegions.PixelIdxList;
        sz = cellfun(@numel,curRegions);
        % first layer filter
        curRegions = curRegions(sz>opts.minSize*opts.minDur/3);
        
        %% identify seed regions
        %% global detect
        seedCandidate = false(numel(curRegions),1);
        t_scl = zeros(numel(curRegions),1);
        fgPix = cell(numel(curRegions),numel(scaleRatios));
        for i = 1:numel(curRegions)
            pix = curRegions{i};
            [ih,iw,it] = ind2sub([H,W,T],pix);
            ihw = unique(sub2ind([H,W],ih,iw));
            dur = max(it)-min(it)+1;
            % second layer filter
            if(dur<opts.minDur || numel(ihw)<opts.minSize)
                continue;
            end 
            labels = setdiff(Map(pix),0);
            if(numel(labels)==1)
                Map(pix) = labels(1);
                continue;
            elseif(numel(labels)>1)
                continue;
            else   % if contain no seed
                seedCandidate(i) = true;
                t_scl(i) = max(1,round(dur/opts.TPatch));
                for j = 1:numel(scaleRatios)
                    pixCur = se.pixResize_half_spa_temp(pix,scaleRatios(j),t_scl(i),H,W,T);
                    fgPix{i,j} = pixCur;
                end
            end
        end

        t_scl = t_scl(seedCandidate);
        seedCandidateRegions = curRegions(seedCandidate);
        fgPix = fgPix(seedCandidate,:);

        sigPass = false(numel(seedCandidateRegions),1);
%         tic;
        parfor i = 1:numel(seedCandidateRegions)
            % - need optimize here -
            [fg,bg1,bg2] = se.neighbor_Individual_spa_temp(fgPix{i,numel(scaleRatios)},fgPix{i,numel(scaleRatios)},t_scl(i),dFResize{end});
            % - need optimize here - 
            % -----------------------------------------
            T_zscore1 = (mean(fg)-mean(bg1))/sqrt(1/numel(fg)+1/numel(bg1));
            T_zscore2 = (mean(fg)-mean(bg2))/sqrt(1/numel(fg)+1/numel(bg2));
            if(T_zscore1<opts.sigThr || T_zscore2<opts.sigThr)
                continue;
            end
            for j = 1:numel(scaleRatios)
                if(~isempty(fgPix{i,j}))
                    % - need optimize here
                    [score1,score2] = se.Individual_Order_Analysis_spa_temp(fgPix{i,j},fgPix{i,j},t_scl(i),dFResize{j});
                    if(score1>opts.sigThr && score2>opts.sigThr)
                        sigPass(i) = true;
                        break;
                    end
                end
            end
        end
%         time = time + toc;

        % Update
        seedCandidateRegions = seedCandidateRegions(sigPass);
        for i = 1:numel(seedCandidateRegions)
            Map(seedCandidateRegions{i}) = nEvt;
            nEvt = nEvt + 1;
        end

    end
    
    sdLst = label2idx(Map);
end