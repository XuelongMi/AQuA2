function evtRun(~,~,f)
% active voxels detection and update overlay map

fh = guidata(f);
opts = getappdata(f,'opts');
if(fh.needTemp.Value && fh.needSpa.Value)
    fprintf('Detecting ...\n')
    try
        opts.sourceSize = str2double(fh.sourceSize.String);
        opts.cDelay = str2double(fh.cDelay.String);
        opts.gtwSmo = str2double(fh.gtwSmo.String);
    %     opts.splitRatio = str2double(fh.splitRatio.String);
        setappdata(f,'opts',opts);
    catch
        msgbox('Error setting parameters')
    end

    seLst1 = getappdata(f,'seLst1');
    subEvtLst1 = getappdata(f,'subEvtLst1');
    dF1 = getappdata(f,'dF1');
    majorInfo1 = getappdata(f,'majorInfo1');
    seLabel1 = getappdata(f,'seLabel1');
    disp('Spatial splitting');
    ff = waitbar(0,'Detecting Channel 1...');
    [subEvtLst1,seLabel1,majorInfo1] = burst.preSplitSub(subEvtLst1,seLabel1,majorInfo1,size(dF1));
    %% pre split by gap
% %     mergingInfo1 = getappdata(f,'mergingInfo1');
% %     ccRegions1 = getappdata(f,'ccRegions1');
% %     dFOrg1 = getappdata(f,'dFOrg1');
% %     [evtLst1,seLabel1] = gapSplit(subEvtLst1,dat1,dF1,majorInfo1,mergingInfo1,ccRegions1,opts);
% %     riseLst1 = [];
% %     datR1 = ones(size(dF1),'uint8')*255;

    [riseLst1,datR1,evtLst1,~] = burst.se2evtTop(dF1,seLst1,subEvtLst1,seLabel1,majorInfo1,opts,ff);
    setappdata(f,'riseLst1',riseLst1);
    setappdata(f,'evt1',evtLst1);

    if(~opts.singleChannel)
        delete(ff);
        dF2 = getappdata(f,'dF2');
        seLst2 = getappdata(f,'seLst2');
        subEvtLst2 = getappdata(f,'subEvtLst2');
        dat2 = getappdata(f,'dat2');
        majorInfo2 = getappdata(f,'majorInfo2');
        seLabel2 = getappdata(f,'seLabel2');
        disp('preprocessing');
        ff = waitbar(0,'Detecting Channel 2...');
        [subEvtLst2,seLabel2,majorInfo2] = preSplitSub(subEvtLst2,seLabel2,majorInfo2,size(dF2));
        [riseLst2,datR2,evtLst2,~] = se2evtTop_alter(dF2,seLst2,subEvtLst2,seLabel2,majorInfo2,opts,ff);
    else
        riseLst2 = []; seLst2 = []; evtLst2 = []; ftsLst2 = []; dffMat2 = []; datR2 = [];
    end
    setappdata(f,'riseLst2',riseLst2);
    setappdata(f,'evt2',evtLst2);
else
    datR1 = [];
    datR2 = [];
    ff = [];
    evtLst1 = getappdata(f,'seLst1');
    evtLst2 = getappdata(f,'seLst2');
    setappdata(f,'evt1',evtLst1);
    setappdata(f,'evt2',evtLst2);
    setappdata(f,'riseLst1',[]);
    setappdata(f,'riseLst2',[]);
end

% setappdata(f,'datR2',datR2);

ui.detect.postRun([],[],f,evtLst1,evtLst2,datR1,datR2,'Events');
% fh.updtFeature1.Enable = 'off';
fh.nEvtName.String = 'nEvt';
if(~opts.singleChannel)
    fh.nEvt.String = [num2str(numel(evtLst1)),' | ',num2str(numel(evtLst2))];
else
    fh.nEvt.String = [num2str(numel(evtLst1))];
end
fprintf('Done\n')
delete(ff);

end





