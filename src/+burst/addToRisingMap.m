function riseLst = addToRisingMap(riseLst,superVoxels,svLabel,dlyMap,nEvt,rgh,rgw,rgt)
% split the rising time map to different events
nEvt0 = max(svLabel(:));
H0 = numel(rgh);
W0 = numel(rgw);
T0 = numel(rgt);
for ii=1:nEvt0
    ids = find(svLabel==ii);
    ihw = [];
    for k = 1:numel(ids)
       pix = superVoxels{ids(k)};
       [ih0,iw0,~] = ind2sub([H0,W0,T0],pix);
        ihw = union(ihw,unique(sub2ind([H0,W0],ih0,iw0)));
    end
    [ihr,iwr] = ind2sub([H0,W0],ihw);
    rghr = min(ihr):max(ihr);
    rgwr = min(iwr):max(iwr);
    mask = false(H0,W0);
    mask(ihw) = true;
    dlyMapr = dlyMap;
    dlyMapr(~mask) = nan;
    dlyMapr = dlyMapr(rghr,rgwr);
    rr = [];
    rr.dlyMap = dlyMapr;
    rr.rgh = min(rgh)+rghr-1;
    rr.rgw = min(rgw)+rgwr-1;
    riseLst{nEvt+ii} = rr;
end
end