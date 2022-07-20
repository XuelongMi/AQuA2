function addCon_tools(f,pTool)
% tools panels
bTool = uix.VBox('Parent',pTool,'Spacing',10);
pLayer = uix.BoxPanel('Parent',bTool,'Title','Layers');
pEvtMngr = uix.BoxPanel('Parent',bTool,'Title','Favourite','Tag','pEvtMngr');
bTool.Heights = [620,-1];

% layer manager
bLayer = uix.VBox('Parent',pLayer,'Padding',3);
ui.com.addConLayer(f,bLayer);

% event manager
bEvt = uix.VBox('Parent',pEvtMngr,'Spacing',3,'Padding',3);
bEvtBtn = uix.HButtonBox('Parent',bEvt,'Spacing',5);
bCurveBtn = uix.HButtonBox('Parent',bEvt,'Spacing',5);
tb = uitable(bEvt,'Data',zeros(0,5),'Tag','evtTable');
bEvtAdd = uix.HBox('Parent',bEvt,'Spacing',10);

uicontrol(bEvtBtn,'String','Select all','Callback',{@ui.evt.evtMngrSelAll,f});
uicontrol(bEvtBtn,'String','Delete','Callback',{@ui.evt.evtMngrDeleteSel,f});
uicontrol(bEvtBtn,'String','Details','Callback',{@ui.evt.showDetails,f});
bEvtBtn.ButtonSize = [120,20];

uicontrol(bCurveBtn,'String','Show curves','Callback',{@ui.evt.evtMngrShowCurve,f});
uicontrol(bCurveBtn,'String','Save curves','Callback',{@ui.evt.saveCurveFig,f});
uicontrol(bCurveBtn,'String','Save waves','Callback',{@ui.evt.saveWaves,f});
% uix.Empty('Parent',bCurveBtn);
bCurveBtn.ButtonSize = [120,20];

tb.ColumnName = {'','CH','Index','Frame','Size','Duration','df/f','Tau'};
tb.ColumnWidth = {15 30 40 40 40 40 40 20};
tb.ColumnEditable = [true,false,false,false,false,false];
tb.CellSelectionCallback = {@ui.evt.evtMngrSelectOne,f};

uicontrol(bEvtAdd,'Style','text','String','Ch1 ID','HorizontalAlignment','left');
uicontrol(bEvtAdd,'Style','edit','Tag','toolsAddEvt1');
uicontrol(bEvtAdd,'String','Add','Callback',{@ui.evt.addOne,f});
uicontrol(bEvtAdd,'Style','text','String','Ch2 ID','HorizontalAlignment','left');
uicontrol(bEvtAdd,'Style','edit','Tag','toolsAddEvt2');
uicontrol(bEvtAdd,'String','Add','Callback',{@ui.evt.addOne2,f});
uix.Empty('Parent',bEvtAdd);
bEvtAdd.Widths = [30,30,40,30,30,40,-1];

bEvt.Heights = [20,20,-1,15];

end








