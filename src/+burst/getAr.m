function [arLst,dActVoxDi] = getAr(dF,opts,evtSpatialMask)
% candidate regions (legacy mode)

% T = size(dF,3);
dActVoxDi = dF>opts.thrARScl;
% for tt=1:T
%     tmp = dF(:,:,tt);
%     tmp = bwareaopen((tmp>opts.thrARScl*sqrt(opts.varEst))&evtSpatialMask,opts.minSize,4);     
% %     tmp = tmp.*evtSpatialMask;
%     dActVoxDi(:,:,tt) = tmp;
% end
arLst = bwconncomp(dActVoxDi);
arLst = arLst.PixelIdxList;

[H,W,T] = size(dF);
sz2D = zeros(numel(arLst),1);
for i = 1:numel(arLst)
    pix = arLst{i};
    [ih,iw,~] = ind2sub([H,W,T],pix);
    ihw = unique(sub2ind([H,W],ih,iw));
    sz2D(i) = numel(ihw);    
end
arLst = arLst(sz2D>=opts.minSize);

end