function [seLstFinal,evtLstFinal,seLabelFInfo,majorityEvtFInfo,opts,sdLst,mergingInfo,ccRegions] = seDetection(dF,dFOrg,arLst,opts,ff)
    [H,W,T] = size(dF);
    activeMap = false(H,W,T);
    for i = 1:numel(arLst)
        activeMap(arLst{i}) = true;
    end
    
    %% setting
    opts.sigma = 1;
    opts.scaleRatios = [7:-1:3];    % useful? necessary????????????? [7:-1:3]
    opts.TPatch = 20;               % naturally it's to avoid too long-duration signal
                                    % downsample the curve to 20 tps
    opts.spaSmo = 3;
    
    tic;
    %% Test Split By my methods from active region
    %% Top-Down
    disp('Resize and smooth');
    dFResize = se.normalizeAndResize(dFOrg,opts);
    
    %% seed detection
    % is there any advantage that when checking seeds, calculate the score
    % for each pixel then combine them together?
    disp('Seed detection');
    [~,Map] = se.seedDetectAc(dF,dFResize,activeMap,opts,ff);
    clear dFResize;
    toc;

    tic;
    %% region grow
    % add temporal penalty for watershed?
    disp('Watershed grow');
    [evtLst,sdLst,ccRegions] = se.markerControlledSplitting_Ac(Map,activeMap,dF,dFOrg,opts);
    if exist('ff','var')
        waitbar(0.7,ff);
    end
    clear Map; clear Map2; 
    toc;
    %% remove empty
    szEvt = cellfun(@numel,evtLst);
    sdLst = sdLst(szEvt>0);
    evtLst = evtLst(szEvt>0);

    %% select propagation part
    tic;
    disp('Majority')
    
    if(~isfield(opts,'major'))
        opts.major = 0.5;
    end
    if(~isfield(opts,'overlap'))
        opts.overlap = 0.5;
    end
    majorityEvt0 = se.getMajority_Ac(sdLst,evtLst,dF);
    toc;
    
    % In temporal dimension, two peaks are closely adjacent, the gap is not obvious?
    if (opts.needRefine)
        disp('Temporally splitting curve according to gap');
        dFSmoVec = reshape(imgaussfilt(dFOrg,4),[],T);
        [evtLst2,majorityEvt2] = se.splitEvt_Ac(dF,dFOrg,evtLst,majorityEvt0,dFSmoVec,opts,1);
        clear dFSmoVec;
    else
        evtLst2 = evtLst;majorityEvt2 = majorityEvt0;
    end

    % Region should further grow?
    %% Grow spatially
    if (opts.needGrow)
        disp('Further grow regions with similar temporal pattern');
        tic;
        [evtLst2,ccRegions] = se.growWatershedResultSpatial(evtLst2,majorityEvt2,dF,opts);
        toc;
    end
    
    %% Merge to super event
    tic;
    disp('Merging signals arising in rising time difference threshold');
    [mergingInfo,majorityEvt2] = se.createMergingInfo(evtLst2,majorityEvt2,ccRegions,dFOrg,opts);
    [seLst,seLstInfoLabel] = se.mergingSEbyInfo(evtLst2,majorityEvt2,mergingInfo,size(dFOrg),ccRegions,opts);
    toc;
    if exist('ff','var')
        waitbar(0.85,ff);
    end
    
    % In temporal dimension, two peaks are closely adjacent, the gap is not obvious?
    %% Merging refine
    if (opts.needRefine)
        tic;
        disp('Merging refine')
        dFSmoVec = reshape(dF,[],T);
        [seLstFinal,evtLstFinal,seLabelFInfo,majorityEvtFInfo,nEvt,nSE,mergingInfo] = se.mergingRefineByInfo(evtLst2,seLstInfoLabel,majorityEvt2,dFOrg,dFSmoVec,ccRegions,mergingInfo,opts);
        toc;   
    else
        seLstFinal = seLst; evtLstFinal = evtLst2; 
        seLabelFInfo = seLstInfoLabel; 
        majorityEvtFInfo = majorityEvt2;
    end
end