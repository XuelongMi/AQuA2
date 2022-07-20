function flow2(~,evtDat,f)
%% the function of RunAllSteps
fh = guidata(f);
opts = getappdata(f,'opts');
setappdata(f,'opts',opts);

ui.detect.preProcessRun([],[],f);
if(~isempty(getappdata(f,'datCorrect1')))
    setappdata(f,'datOrg1',getappdata(f,'datCorrect1'));
    rmappdata(f,'datCorrect1');
    setappdata(f,'datOrg2',getappdata(f,'datCorrect2'));
    rmappdata(f,'datCorrect2');
end
fh.registrateCorrect.Enable = 'off';
fh.bleachCorrect.Enable = 'off';

ui.detect.actRun([],[],f);

ui.detect.phaseRun([],[],f);

ui.detect.evtRun([],[],f);

ui.detect.gloRun([],[],f);

ui.detect.feaRun([],[],f);

% controls
fh.deOutBack.Visible = 'on';
fh.deOutRun.String = 'Extract';
fh.deOutNext.Visible = 'off';

fh.deOutTab.Selection = 6;
for i=1:6
    fh.deOutTab.TabEnables{i} = 'on';
end

end


