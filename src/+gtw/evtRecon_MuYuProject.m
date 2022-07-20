function [evtL,evtRecon] = evtRecon_MuYuProject(spLst,cx,evtMap0,minShow,seMapValid)

[H,W] = size(evtMap0);
[nSp,T] = size(cx);

seMapValid = reshape(seMapValid,[],T);
evtRecon = zeros(H*W,T);
evtL = zeros(H*W,T);
% seedMap1 = zeros(H,W);
for ii=1:nSp
    sp0 = spLst{ii};
    x0 = cx(ii,:);
    evtRecon(sp0,:) = repmat(x0,numel(sp0),1);
    l0 = mode(evtMap0(sp0));
    t0 = find(x0>=minShow,1);
    t1 = find(x0>=minShow,1,'last');
    
    spLabels = zeros(numel(sp0),T);
    spLabels(:,t0:t1) = l0;
    spLabels(seMapValid(sp0,:)) = l0;
    evtL(sp0,:) = spLabels;
%     [ih,iw] = ind2sub([H,W],sp0);
%     seedMap1(round(mean(ih)),round(mean(iw))) = ii;
end
evtRecon = reshape(evtRecon,H,W,T);
evtL = reshape(evtL,H,W,T);

end
