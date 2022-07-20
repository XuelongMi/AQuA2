function flow(~,evtDat,f,op)

fh = guidata(f);
nTabTot = numel(fh.deOutTab.TabEnables);
ixTab = fh.deOutTab.Selection;
opts = getappdata(f,'opts');
setappdata(f,'opts',opts);

% controls
if strcmp(op,'chg')
    ixTab = evtDat.NewValue;
    fh.deOutBack.Visible = 'on';
    fh.deOutNext.Visible = 'on';
    fh.deOutNext.Enable = 'on';
    fh.deOutRun.String = 'Run';
    fh.deOutNext.String = 'Next';
    switch ixTab
        case 1
            fh.deOutBack.Visible = 'off';
        case 6
            fh.deOutRun.String = 'Extract';
%             fh.deOutNext.Visible = 'off';
            fh.deOutNext.String = 'CFU detect';
            if(isempty(getappdata(f,'finishAll')))
                fh.deOutNext.Enable = 'off';
            end
%             fh.deOutNext.Visible = 'off';
    end
    if ixTab<nTabTot
        if strcmp(fh.deOutTab.TabEnables{ixTab},'off')
            fh.deOutNext.Enable = 'off';
        end
    end
end

% go to previous step
if strcmp(op,'back')
    if ixTab>1
        fh.deOutTab.Selection = ixTab-1;
    end
    fh.deOutNext.Enable = 'on';
end

% run current step
if strcmp(op,'run')
    switch ixTab
        case 1
            ui.detect.preProcessRun([],[],f);
        case 2
            ui.detect.actRun([],[],f);
        case 3
            ui.detect.phaseRun([],[],f);
        case 4
            ui.detect.evtRun([],[],f);
        case 5
            ui.detect.gloRun([],[],f);
        case 6
            ui.detect.feaRun([],[],f);
    end
    if ixTab<nTabTot
        fh.deOutNext.Enable = 'on';
    end
end

% go to next step
if strcmp(op,'next')
    if ixTab<nTabTot
        if(fh.deOutTab.Selection==1 && strcmp(fh.registrateCorrect.Enable,'on'))
            selection = questdlg('Use current processed data? Enter next step, the preprocessing setting cannot be changed?', ...
                    'warning','OK','Cancel','Cancel');
            switch selection
                case 'OK'
                    if(~isempty(getappdata(f,'datCorrect1')))
                        setappdata(f,'datOrg1',getappdata(f,'datCorrect1'));
                        rmappdata(f,'datCorrect1');
                        setappdata(f,'datOrg2',getappdata(f,'datCorrect2'));
                        rmappdata(f,'datCorrect2');
                    end
                    fh.registrateCorrect.Enable = 'off';
                    fh.bleachCorrect.Enable = 'off';
                    
                case 'Cancel'
                    return
            end
        end
        
        fh.deOutTab.Selection = ixTab+1;
        fh.deOutTab.TabEnables{ixTab+1} = 'on';

        
        
        if(fh.deOutTab.Selection==5)
            fh.mergeEventCorr.Visible = 'on';
            fh.mergeEventCorrText.Visible = 'on';
        end
    else
        ui.detect.CFURun([],[],f);
    end
end
end


