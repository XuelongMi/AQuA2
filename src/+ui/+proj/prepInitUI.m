function prepInitUI(f,fh,opts,scl,~,stg,~)
    
    % layer panel
    gap0 = [0.01 0.1];
    fh.sldMin.Min = scl.min;
    fh.sldMin.Max = scl.max;
    fh.sldMin.SliderStep = gap0;
    fh.sldMin.Value = scl.min;
    
    fh.sldMax.Min = scl.min;
    fh.sldMax.Max = scl.max;
    fh.sldMax.SliderStep = gap0;
    fh.sldMax.Value = scl.max;
    
    fh.sldBri1.Min = 0.1;
    fh.sldBri1.Max = 10;
    fh.sldBri1.SliderStep = gap0;
    fh.sldBri1.Value = scl.bri1;
    
    fh.sldBri2.Min = 0.1;
    fh.sldBri2.Max = 10;
    fh.sldBri2.SliderStep = gap0;
    fh.sldBri2.Value = scl.bri1;
    
    fh.sldBriL.Min = 0.1;
    fh.sldBriL.Max = 10;
    fh.sldBriL.SliderStep = gap0;
    fh.sldBriL.Value = scl.briL;
    
    fh.sldBriR.Min = 0.1;
    fh.sldBriR.Max = 10;
    fh.sldBriR.SliderStep = gap0;
    fh.sldBriR.Value = scl.briR;
    
    fh.sldMinOv.Min = 0;
    fh.sldMinOv.Max = 1;
    fh.sldMinOv.SliderStep = gap0;
    fh.sldMinOv.Value = scl.minOv;
    
    fh.sldMaxOv.Min = 0;
    fh.sldMaxOv.Max = 1;
    fh.sldMaxOv.SliderStep = gap0;
    fh.sldMaxOv.Value = scl.maxOv;
    
    fh.sldBriOv.Min = 0;
    fh.sldBriOv.Max = 1;
    fh.sldBriOv.SliderStep = gap0;
    fh.sldBriOv.Value = scl.briOv;
    
    % data panel
    fh.sldMov.Min = 1;
    fh.sldMov.Max = opts.sz(3);
    fh.sldMov.SliderStep = [1/opts.sz(3),0.05];
%     fh.sldMov.BlockIncrement = 1;
%     fh.sldMov.VisibleAmount = 0;
    fh.sldMov.Value = 1;
    
    % detection parameters
    try
        fh.bleachCorrect.Value = opts.bleachCorrect;
        if(fh.bleachCorrect.Value==0)
            fh.bleachCorrect.Value = 1;
        end
    catch
        fh.bleachCorrect.Value = 1;
    end
    
    fh.thrArScl.String = num2str(opts.thrARScl);
    fh.smoXY.String = num2str(opts.smoXY);
    fh.minSize.String = num2str(opts.minSize);
    try
        fh.smoT.String = num2str(opts.smoT);
    catch
        
    end
    
    fh.minDur.String = num2str(opts.minDur);
    fh.ratio.String = num2str(opts.ratio);
    fh.sigThr.String = num2str(opts.sigThr);
    try
        fh.maxDelay.String = num2str(opts.maxDelay);
    catch
        
    end
    
    try
        fh.splitRatio.String = num2str(opts.splitRatio);
    catch
        
    end
    fh.cRise.String = num2str(opts.cRise);
    fh.cDelay.String = num2str(opts.cDelay);
    fh.gtwSmo.String = num2str(opts.gtwSmo);
    
    fh.zThr.String = num2str(opts.zThr);
    
    fh.ignoreMerge.Value = 1*(opts.ignoreMerge>0);
    fh.mergeEventDiscon.String = num2str(opts.mergeEventDiscon);
    fh.mergeEventCorr.String = num2str(opts.mergeEventCorr);
    fh.mergeEventMaxTimeDif.String = num2str(opts.mergeEventMaxTimeDif);
    
    fh.extendEvtRe.Value = 1*(opts.extendEvtRe>0);
    
    fh.ignoreTau.Value = 1*(opts.ignoreTau>0);
    
    % color overlay
    ui.over.getColMap([],[],f);
    
    try
        % update overlay menu
        ui.over.updateOvFtMenu([],[],f);
        
        % User defined features
        ui.over.chgOv([],[],f,0);
        ui.over.chgOv([],[],f,1);
        ui.over.chgOv([],[],f,2);
        ui.evt.evtMngrRefresh([],[],f);
    catch
    end
    
    % resize GUI
    fh.g.Selection = 3;
    f.Resize = 'on';
    f.Position = getappdata(f,'guiMainSz');
    
    dbgx = getappdata(f,'dbg');
    if isempty(dbgx); dbgx=0; end        
    
    % UI visibility according to steps
    if stg.detect==0  % not started yet
        xx = fh.deOutTab.TabEnables;
        for ii=2:numel(xx)
            xx{ii} = 'off';
        end
        fh.deOutTab.TabEnables = xx;
        fh.deOutNext.Enable = 'off';
        fh.pFilter.Visible = 'off';
        fh.pExport.Visible = 'off';
        fh.pEvtMngr.Visible = 'off';
        fh.pSys.Visible = 'off';
        fh.deOutTab.Selection = 1;
        fh.deOutBack.Visible = 'off';
    else  % finished
        ui.detect.filterInit([],[],f);
        evtLst = getappdata(f,'evtLst');
        xx = fh.deOutTab.TabEnables;
        for ii=1:numel(xx)-1
            xx{ii} = 'off';
        end
        fh.deOutTab.TabEnables = xx;
        fh.deOutBack.Enable = 'off';
        fh.deOutTab.Selection = numel(xx);
%         if isempty(evtLst) && dbgx==0
%             fh.bWkfl.Heights(2) = 0;  % never show detection part again
%         end
    end

    
    % show movie
    ui.movStep(f,1,[],1);
    
end





