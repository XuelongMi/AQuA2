function actRun(~,~,f)
% active voxels detection and update overlay map

fprintf('Detecting ...\n')

fh = guidata(f);
bd = getappdata(f,'bd');
opts = getappdata(f,'opts');

% only inside user drawn cells
sz = opts.sz;
evtSpatialMask = true(sz(1),sz(2));
if bd.isKey('cell')
    bd0 = bd('cell');
    evtSpatialMask = false(sz(1),sz(2));
    for ii=1:numel(bd0)
        p0 = bd0{ii}{2};
        evtSpatialMask(p0) = true;
    end
end

ff = waitbar(0,'Removing baseline + Noise estimation...');
if(~isfield(opts,'alreadyRemoveBaseline') || ~opts.alreadyRemoveBaseline ||...
        opts.noiseEstimation~=fh.noiseEstimation.Value || opts.smoXY~= str2double(fh.smoXY.String))
    opts.noiseEstimation = fh.noiseEstimation.Value;
    opts.smoXY = str2double(fh.smoXY.String);
    xx = randn(10000,opts.cut); %% bias need to be reconsidered
    xxMA = movmean(xx,opts.movAvgWin,2);
    xxMin = min(xxMA,[],2);
    opts.xBias = nanmean(xxMin(:));
    datOrg1 = getappdata(f,'datOrg1');
    [dat1,dF1,dFOrg1,opts] = burst.baselineRemoveAndNoiseEstimation(datOrg1,opts,evtSpatialMask,ff);
    opts.maxdF1 = quantile(dF1(:),0.999);
    setappdata(f,'dat1',dat1);
    setappdata(f,'dF1',dF1);
    setappdata(f,'dFOrg1',dFOrg1);
    clear dat1 dF1 dFOrg1 datOrg1;
    waitbar(0.25,ff);
    if(~opts.singleChannel)
        datOrg2 = getappdata(f,'datOrg2');
        [dat2,dF2,dFOrg2,opts] = burst.baselineRemoveAndNoiseEstimation(datOrg2,opts,evtSpatialMask,ff);
        opts.maxdF2 = quantile(dF2(:),0.999);
    else
        dat2 = []; dF2 = []; dFOrg2 = []; arLst2 = [];
    end
    setappdata(f,'dat2',dat2);
    setappdata(f,'dF2',dF2);
    setappdata(f,'dFOrg2',dFOrg2);
    
    clear dat2 dF2 dFOrg2 datOrg2;
    opts.alreadyRemoveBaseline = true;
    setappdata(f,'opts',opts);
end
waitbar(0.5,ff,'Detecting active regions...');

opts.thrARScl = str2double(fh.thrArScl.String);
opts.minSize = str2double(fh.minSize.String);
opts.maxSize = str2double(fh.maxSize.String);
opts.minDur = str2double(fh.minDur.String);
opts.circularityThr = str2double(fh.circularityThr.String);
  
opts.compress = 0.3;
dF1 = getappdata(f,'dF1');
[arLst1] = burst.acDetect(dF1,opts,evtSpatialMask,ff);  % foreground and seed detection
setappdata(f,'arLst1',arLst1);
clear dF1;
if(~opts.singleChannel)
    dF2 = getappdata(f,'dF2');
    [arLst2] = burst.acDetect(dF2,opts,evtSpatialMask,ff);  % foreground and seed detection
    clear dF2;
    
else
    arLst2 = [];
end
setappdata(f,'arLst2',arLst2);

waitbar(1,ff);

ui.detect.postRun([],[],f,arLst1,arLst2,[],[],'Step 1: active voxels');
fh.GaussFilter.Enable = 'on';


fh.nEvtName.String = 'nActiveRegions';
if(~opts.singleChannel)
    fh.nEvt.String = [num2str(numel(arLst1)),' | ',num2str(numel(arLst2))];
else
    fh.nEvt.String = [num2str(numel(arLst1))];
end
fprintf('Done\n');
delete(ff);
end