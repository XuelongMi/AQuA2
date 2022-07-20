function loadOpt(~,~,f)
    
    %file0 = [opts.fileName,'_AQuA']; SP, 18.07.16
    opts = getappdata(f,'opts');
    [file,path] = uigetfile('.csv','Choose Parameter file',opts.filePath1);
    if ~isnumeric([path,file])
        optsOrg = getappdata(f,'opts');
        opts = ui.proj.csv2struct([path,file]);
        opts.filePath1 = optsOrg.filePath1;
        opts.filePath2 = optsOrg.filePath2;
        opts.fileName1 = optsOrg.fileName1;
        opts.fileName2 = optsOrg.fileName2;
        opts.fileNameType1 = optsOrg.fileType1;
        opts.fileNameType2 = optsOrg.fileType2;
        opts.sz = optsOrg.sz;
        opts.maxValueDat = optsOrg.maxValueDat;
        opts.maxValueDepth = optsOrg.maxValueDepth;
        opts.frameRate = optsOrg.frameRate;
        opts.spatialRes = optsOrg.spatialRes;
        opts.regMaskGap = optsOrg.regMaskGap;
        setappdata(f,'opts',opts);
        
        
        % adjust interface parameters
        fh = guidata(f);
        fh.bleachCorrect.Value = opts.bleachCorrect;
        
        fh.thrArScl.String = num2str(opts.thrARScl);
        fh.smoXY.String = num2str(opts.smoXY);
        fh.minSize.String = num2str(opts.minSize);
        fh.smoT.String = num2str(opts.smoT);
        
        fh.minDur.String = num2str(opts.minDur);
        fh.ratio.String = num2str(opts.ratio);
        fh.sigThr.String = num2str(opts.sigThr);
        fh.maxDelay.String = num2str(opts.maxDelay);
        
        fh.cRise.String = num2str(opts.cRise);
        fh.cDelay.String = num2str(opts.cDelay);
        fh.gtwSmo.String = num2str(opts.gtwSmo);
        
        fh.zThr.String = num2str(opts.zThr);
        
        fh.ignoreMerge.Value = opts.ignoreMerge;
        fh.mergeEventDiscon.String = num2str(opts.mergeEventDiscon);
        fh.mergeEventCorr.String = num2str(opts.mergeEventCorr);
        fh.mergeEventMaxTimeDif.String = num2str(opts.mergeEventMaxTimeDif);
        
        fh.extendEvtRe.Value = opts.extendEvtRe;
        
        fh.ignoreTau.Value = opts.ignoreTau;
        
        n = round(fh.sldMov.Value);
        ui.mov.updtMovInfo(f,n,opts.sz(3));
    end
end