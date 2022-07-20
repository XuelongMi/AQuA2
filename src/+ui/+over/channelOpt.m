function channelOpt(~,~,f)
    btSt = getappdata(f,'btSt');
    fh = guidata(f);
    btSt.channelSelect = fh.channelOption.Value;
    setappdata(f,'btSt',btSt);
    ui.movStep(f);

end