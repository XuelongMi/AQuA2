function addDetectTab(f,pDeOut)
    % addDetectTab adds event detection pipeline panels in tabs    
    h = cell(0);
    
    % top
    bDeOut = uix.VBox('Parent',pDeOut,'Padding',3,'Spacing',3);
    uix.Empty('Parent',bDeOut);
    deOutTab = uix.TabPanel('Parent',bDeOut,'Tag','deOutTab');
    deOutCon = uix.HButtonBox('Parent',bDeOut,'Spacing',10);
    deOutRunAll = uix.HButtonBox('Parent',bDeOut,'Spacing',10);
    uix.Empty('Parent',bDeOut);
    bDeOut.Heights = [-1,260,25,25,-1];
    
    % tabs
    pPre = uix.Panel('Parent',deOutTab,'Title','Preprocess','Tag','pBc');
    pAct = uix.Panel('Parent',deOutTab,'Title','Active region detection','Tag','pAct');
    pSp = uix.Panel('Parent',deOutTab,'Title','Temporal segmentation => super events','Tag','pSp');
    pEvt = uix.Panel('Parent',deOutTab,'Title','Spatial segmentation => events','Tag','pEvt');
    pGlo = uix.Panel('Parent',deOutTab,'Title','Global signal detection','Tag','pGlo');
    pFea = uix.Panel('Parent',deOutTab,'Title','Feature extraction','Tag','pFea');
    deOutTab.TabTitles = {'Preprocess','Active region','Temporal','Spatial','Global','Feature'};
    deOutTab.TabWidth = 40;
    deOutTab.SelectionChangedFcn = {@ui.detect.flow,f,'chg'};
    
    % controls
    pBack = uicontrol(deOutCon,'String','Back','Tag','deOutBack');
    pRun = uicontrol(deOutCon,'String','Run','Tag','deOutRun');
    pNext = uicontrol(deOutCon,'String','Next','Tag','deOutNext');
    pBack.Callback = {@ui.detect.flow,f,'back'};
    pRun.Callback = {@ui.detect.flow,f,'run'};
    pNext.Callback = {@ui.detect.flow,f,'next'};
    deOutCon.ButtonSize = [100,20];
    pSaveOpt = uicontrol(deOutRunAll,'String','SaveOpts','Tag','deSaveOpt');
    pLoadOpt = uicontrol(deOutRunAll,'String','LoadOpts','Tag','deLoadOpt');
    pRunAll = uicontrol(deOutRunAll,'String','RunAllSteps','Tag','deOutRunAll');
    pSaveOpt.Callback = {@ui.detect.saveOpt,f};
    pSaveOpt.BackgroundColor = [.3,.5,.8];
    pSaveOpt.ForegroundColor = [1,1,1];
    pLoadOpt.Callback = {@ui.detect.loadOpt,f};
    pLoadOpt.BackgroundColor = [.3,.5,.8];
    pLoadOpt.ForegroundColor = [1,1,1];
    pRunAll.Callback = {@ui.detect.flow2,f};
    pRunAll.BackgroundColor = [.3,.5,.8];
    pRunAll.ForegroundColor = [1,1,1];
    deOutRunAll.ButtonSize = [100,20];
    
    % Bleach Correction
    bPre = uix.VBox('Parent',pPre);
    uicontrol(bPre,'Style','text','String','--------------Registration--------------');
    uicontrol(bPre,'Style','popupmenu','Value',1,'Tag','registrateCorrect','String',{'Not registrate','Rigid registration by cross correlation based on channel 1','Rigid registration by cross correlation based on channel 2'});
    uicontrol(bPre,'Style','text','String','--------------Photobleach correction--------------');
    uicontrol(bPre,'Style','popupmenu','Value',1,'Tag','bleachCorrect','String',{'Not correct bleach','Remove bleach globally','Remove bleach by intensity'});
    uix.Empty('Parent',bPre);
    bPre.Heights = [20,20,20,20,-1];
    
    % event detection: active voxels
    bAct = uix.VBox('Parent',pAct);
    uicontrol(bAct,'Style','text','String','--------------Smoothing--------------');
    gAct1 = uix.Grid('Parent',bAct,'Padding',0,'Spacing',8);
    uicontrol(gAct1,'Style','edit','String','0.5','Tag','smoXY');
    uicontrol(gAct1,'Style','text','String','Spatial smoothing (sigma)');
    gAct1.Widths = [50,-1]; gAct1.Heights = [15];
    uicontrol(bAct,'Style','text','String','--------------Noise estimation--------------');
    uicontrol(bAct,'Style','popupmenu','Value',1,'Tag','noiseEstimation','String',{'Uniform noise estimation for whole video','Pixel-wise noise estimation'});
    uicontrol(bAct,'Style','text','String','--------------Region detection setting--------------');
    gAct2 = uix.Grid('Parent',bAct,'Padding',0,'Spacing',8);
    uicontrol(gAct2,'Style','edit','String','3','Tag','thrArScl');
    uicontrol(gAct2,'Style','edit','String','5','Tag','minDur');
    uicontrol(gAct2,'Style','edit','String','8','Tag','minSize');
    uicontrol(gAct2,'Style','edit','String','10000','Tag','maxSize');
    uicontrol(gAct2,'Style','edit','String','0','Tag','circularityThr');
    % -----------------------------------------------------------------------
    uicontrol(gAct2,'Style','text','String','Intensity threshold scaling factor');
    uicontrol(gAct2,'Style','text','String','Minimum duration');
    uicontrol(gAct2,'Style','text','String','Minimum size (pixels)');
    uicontrol(gAct2,'Style','text','String','Maximum size (pixels)');
    uicontrol(gAct2,'Style','text','String','Circurlarity threshold for active region');
    gAct2.Widths = [50,-1]; gAct2.Heights = [15,15,15,15,15];
    bAct.Heights = [25,20,25,20,25,-1];
    
    % event detection: superpixels and rising time
    bPhase = uix.VBox('Parent',pSp);
    uicontrol(bPhase,'Style','checkbox','String',{'Need temporal segmentation?'},...
        'Value',1,'Tag','needTemp','Callback',{@ui.com.ableTemp,f});
    uicontrol(bPhase,'Style','text','String','--------------Temporal segmentation setting--------------','Tag','tempString');
    gPhase = uix.Grid('Parent',bPhase,'Padding',10,'Spacing',8,'Tag','tempSetting');
    uicontrol(gPhase,'Style','edit','String','3.5','Tag','sigThr');
    uicontrol(gPhase,'Style','edit','String','0.8','Tag','maxDelay');
    h{end+1} = uicontrol(gPhase,'Style','text','String','Zscore of seed detection');
    h{end+1} = uicontrol(gPhase,'Style','text','String','Allowed signal difference in SE');
    uicontrol(bPhase,'Style','checkbox','String',{'No obvious gap in two temporally adjacent signals?'},...
        'Value',0,'Tag','needRefine');
    uicontrol(bPhase,'Style','checkbox','String',{'Are active regions not large enough?'},...
        'Value',0,'Tag','needGrow');
    gPhase.Widths = [50,-1]; gPhase.Heights = [15,15];
    bPhase.Heights = [20,20,50,20,20];
    
    % event detection: events
    bEvt = uix.VBox('Parent',pEvt);
    uicontrol(bEvt,'Style','checkbox','String',{'Need spatial segmentation?'},...
        'Value',1,'Tag','needSpa','Callback',{@ui.com.ableSpa,f});
    uicontrol(bEvt,'Style','text','String','--------------Spatial segmentation setting--------------','Tag','spaString');
    gEvt = uix.Grid('Parent',bEvt,'Padding',10,'Spacing',8,'Tag','spaSetting');
    uicontrol(gEvt,'Style','edit','String','200','Tag','sourceSize');
    uicontrol(gEvt,'Style','edit','String','0.5','Tag','cDelay');
    uicontrol(gEvt,'Style','edit','String','0.1','Tag','gtwSmo');
    
    h{end+1} = uicontrol(gEvt,'Style','text','String','Minimum source size');
    h{end+1} = uicontrol(gEvt,'Style','text','String','Slowest delay in propagation');
    h{end+1} = uicontrol(gEvt,'Style','text','String','Propagation smoothness');
    gEvt.Widths = [50,-1]; gEvt.Heights = [15,15,15];
    bEvt.Heights = [20,20,80];
    
    % Global signal
    bGlo = uix.VBox('Parent',pGlo,'Padding',10);
    uicontrol(bGlo,'Style','checkbox','String','Detect global signal',...
        'Value',0,'Tag','detectGlo');
    uix.Empty('Parent',bGlo);
    bGlo.Heights = [20,-1];
    
    % Extract feature
    bFea = uix.VBox('Parent',pFea,'Padding',10);
    uicontrol(bFea,'Style','checkbox','String','Ignore delay Tau',...
        'Value',1,'Tag','ignoreTau');
    uix.Empty('Parent',bFea);
    bFea.Heights = [20,-1];
    
    for ii=1:numel(h)
        h00 = h{ii};
        h00.HorizontalAlignment = 'left';
    end    
end

