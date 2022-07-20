function majorityEvt0 = getMajority_Ac(sdLst,evtLst,dF)
    [H,W,T] = size(dF);
    majorityEvt0 = cell(numel(sdLst),1);
    dF = reshape(dF,[],T);
    parfor i = 1:numel(sdLst)
        % initialization
        [ih,iw,it] = ind2sub([H,W,T],sdLst{i});
        t0 = min(it);
        t1 = max(it);
        ihw = unique(sub2ind([H,W],ih,iw));
        % seed to majority
        [mIhw,TW,delays] = se.seed2Majoirty(ihw,dF,[H,W,T],evtLst{i},t0:t1);
        % update
        majorityEvt0{i}.ihw = mIhw;
        majorityEvt0{i}.TW = TW;
        majorityEvt0{i}.delays = delays;
        majorityEvt0{i}.needUpdatePeak = true;
        
        curve00 = mean(dF(mIhw,TW),1);
        curve00 = curve00 - min(curve00);
        curve00 = curve00/max(curve00);
        majorityEvt0{i}.curve = curve00;
        
    end
end