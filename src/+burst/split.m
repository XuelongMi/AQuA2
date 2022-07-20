function lblMapS = split(lblMap,dF)
%% some super voxels have multiple peak. 
    svLst = label2idx(lblMap);
    [H,W,T] = size(lblMap);
    dFVec = reshape(dF,[],T);
    lblMapS = zeros(H,W,T);
    cnt = 1;
    for i = 1:numel(svLst)
        pix = svLst{i};
        [ih,iw,it] = ind2sub([H,W,T],pix);
        rgT = min(it):max(it);
        ihw = sub2ind([H,W],ih,iw);
        curve = mean(dFVec(ihw,:));
        s1 = sqrt(median((curve(2:end)-curve(1:end-1)).^2)/0.9133);
        TW = burst.getPeak(curve,rgT,s1);
        
        if isempty(TW)
           TW = cell(1);
           TW{1} = rgT;
        end
        
        for j = 1:numel(TW)
            rgT0 = TW{j};
            it0 = min(rgT0);
            it1 = max(rgT0);
            l = find(it>=it0 & it<=it1);
            lblMapS(pix(l)) = cnt;
            cnt = cnt+1;
        end
    end

end