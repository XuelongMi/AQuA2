function ftsBase = getBasicFeatures(voli0,muPerPix,nEvt,ftsBase)
% getFeatures extract local features from events

% basic features
ftsBase.map{nEvt} = sum(voli0,3);
cc = regionprops(ftsBase.map{nEvt}>0,'Area');
ftsBase.area(nEvt) = sum([cc.Area])*muPerPix*muPerPix;

% set 50% to remove artifact influence
cc = regionprops(ftsBase.map{nEvt}>0.5*max(ftsBase.map{nEvt}(:)),'Area','Perimeter');
ftsBase.peri(nEvt) = sum([cc.Perimeter]);
ftsBase.circMetric(nEvt) = (sum([cc.Perimeter])^2/(4*pi*sum([cc.Area])));

end









