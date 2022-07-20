function [seLst,evtLst,seLabel,majorityEvt,nEvt,nSE,mergingInfo] = mergingRefineByInfo(evtLst,seLabel,majorityEvt,dFOrg,dFSmoVec,ccRegions,mergingInfo,opts)

    MaxIter = 10;
    nEvt = zeros(MaxIter+1,1);
    nSE = zeros(MaxIter+1,1);
    if(numel(evtLst)==0)
        seLst = cell(1,0);
        evtLst = cell(1,0);
        seLabel = [];
        majorityEvt = cell(1,0);
        nEvt = 0;
        nSE = 0;
        return;
    end
    nEvt(1) = numel(evtLst);
    nSE(1) = max(seLabel);
    
    seLst = cell(max(seLabel),1);
    for i = 1:numel(evtLst)
       seID = seLabel(i);
       seLst{seID} = [seLst{seID};evtLst{i}];
    end    
    

    for k = 1:MaxIter
        [evtLst,majorityEvt,mergingInfo,refineWork] = se.refineEvtsByInfo2(evtLst,seLabel,majorityEvt,dFOrg,dFSmoVec,mergingInfo,opts);
        if(~refineWork)
            break;
        end
        mergingInfo = se.updateMergingInfo(evtLst,dFOrg,majorityEvt,mergingInfo);
        [seLst,seLabel,mergingInfo] = se.mergingSEbyInfo(evtLst,majorityEvt,mergingInfo,size(dFOrg),ccRegions,opts);
        nEvt(k+1) = numel(evtLst);
        nSE(k+1) = numel(seLst);
    end
    
    [evtLst,majorityEvt,mergingInfo,refineWork] = se.refineEvtsCorrectByInfo(evtLst,seLabel,majorityEvt,dFOrg,dFSmoVec,mergingInfo,opts);
    if(refineWork)
        mergingInfo = se.updateMergingInfo(evtLst,dFOrg,majorityEvt,mergingInfo);
        seLabel = [seLabel,(numel(seLabel)+1):numel(evtLst)];
        [evtLst,majorityEvt,mergingInfo,refineWork] = se.refineEvtsByInfo2(evtLst,seLabel,majorityEvt,dFOrg,dFSmoVec,mergingInfo,opts);
        mergingInfo = se.updateMergingInfo(evtLst,dFOrg,majorityEvt,mergingInfo);
        [seLst,seLabel,mergingInfo] = se.mergingSEbyInfo(evtLst,majorityEvt,mergingInfo,size(dFOrg),ccRegions,opts);
    end
end