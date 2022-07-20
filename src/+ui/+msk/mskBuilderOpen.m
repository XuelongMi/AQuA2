function mskBuilderOpen(~,~,f)
    fh = guidata(f);
    
    fgOK = 0;
    bgOK = 0;
    opts = getappdata(f,'opts');
    bd = getappdata(f,'bd');
    if bd.isKey('maskLst')
        bdMsk = bd('maskLst');
        for ii=1:numel(bdMsk)
            rr = bdMsk{ii};
            if strcmp(rr.type,'Channel1')
                fgOK = 1;
            end
            if strcmp(rr.type,'Channel2')
                bgOK = 1;
            end
        end
    end
    
    if fgOK==0
        ui.msk.readMsk([],[],f,'self_CH1','Channel1',0);
    end
    if (bgOK==0 & ~opts.singleChannel)
        ui.msk.readMsk([],[],f,'self_CH2','Channel2',0);
    end
    
    ui.msk.mskLstViewer([],[],f,'refresh');
    
    fh.g.Selection = 4;
end