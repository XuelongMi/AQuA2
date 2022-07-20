function stepOne(~,~,f)
fh = guidata(f);
n = round(fh.sldMov.Value);
ui.movStep(f,n,[]);
if(isfield(fh,'showcurves') && ~isempty(fh.showcurves))
    evtIdx = fh.showcurves(:,1);
    if ~isempty(evtIdx)
        channels = fh.showcurves(:,2);
        evtIdx1 = evtIdx(channels==1);
        evtIdx2 = evtIdx(channels==2);
        ui.evt.curveRefresh([],[],f,evtIdx1,evtIdx2);
    end
end
end