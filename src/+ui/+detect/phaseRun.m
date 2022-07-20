function phaseRun(~,~,f)
% active voxels detection and update overlay map

fprintf('Detecting ...\n')

fh = guidata(f);

if(fh.needTemp.Value)
    opts = getappdata(f,'opts');
    opts.needTemp = fh.needTemp.Value;
    opts.ratio = 0.5;
    opts.sigThr = str2double(fh.sigThr.String);
    opts.maxDelay = str2double(fh.maxDelay.String);
    opts.needRefine = fh.needRefine.Value;
    opts.needGrow = fh.needGrow.Value;

    dF1 = getappdata(f,'dF1');
    dFOrg1 = getappdata(f,'dFOrg1');
    arLst1 = getappdata(f,'arLst1');
    ff = waitbar(0,'Detecting Channel1 ...');
    [seLst1,subEvtLst1,seLabel1,majorInfo1,opts,sdLst1,~,~] = se.seDetection(dF1,dFOrg1,arLst1,opts,ff);

    if(~opts.singleChannel)
        delete(ff);
        dF2 = getappdata(f,'dF2');
        dFOrg2 = getappdata(f,'dFOrg2');
        arLst2 = getappdata(f,'arLst2');
        ff = waitbar(0,'Detecting Channel2 ...');
        [seLst2,subEvtLst2,seLabel2,majorInfo2,opts,sdLst2,~,~] = se.seDetection(dF2,dFOrg2,arLst2,opts,ff);
    else
        seLst2 = [];
        sdLst2 = [];
        subEvtLst2 = [];
        seLabel2 = [];
        majorInfo2 = [];
    end

    % save data
    setappdata(f,'subEvtLst1',subEvtLst1);
    setappdata(f,'seLst1',seLst1);
    setappdata(f,'seLabel1',seLabel1);
    setappdata(f,'majorInfo1',majorInfo1);
    setappdata(f,'subEvtLst2',subEvtLst2);
    setappdata(f,'seLst2',seLst2);
    setappdata(f,'seLabel2',seLabel2);
    setappdata(f,'majorInfo2',majorInfo2);
    setappdata(f,'opts',opts);
else
    opts = getappdata(f,'opts');
    arLst1 = getappdata(f,'arLst1');
    if(~opts.singleChannel)
        arLst2 = getappdata(f,'arLst2');
    else
        arLst2 = [];
    end
    sdLst1 = arLst1;
    sdLst2 = arLst2;
    subEvtLst1 = arLst1;
    subEvtLst2 = arLst2;
    seLst1 = arLst1;
    seLst2 = arLst2;
    setappdata(f,'seLst1',seLst1);
    setappdata(f,'seLst2',seLst2);
    ff = [];
end

ui.detect.postRun([],[],f,sdLst1,sdLst2,[],[],'Step 2aa: seeds');
ui.detect.postRun([],[],f,subEvtLst1,subEvtLst2,[],[],'Step 2a: watershed results');
ui.detect.postRun([],[],f,seLst1,seLst2,[],[],'Step 2b: super events');

fh.nEvtName.String = 'nSe';
if(~opts.singleChannel)
    fh.nEvt.String = [num2str(numel(seLst1)),' | ',num2str(numel(seLst2))];
else
    fh.nEvt.String = [num2str(numel(seLst1))];
end

fprintf('Done\n')
delete(ff);
end





