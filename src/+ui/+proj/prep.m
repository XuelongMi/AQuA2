function prep(~,~,f,op,res)
% read data or load experiment
% op
% 0: new project
% 1: load project or saved results
% 2: load from workspace
% FIXME: udpate GUI settings (btSt), instead of re-build it

fprintf('Loading ...\n');
ff = waitbar(0,'Loading ...');

% cfgFile = 'uicfg.mat';
% if ~exist(cfgFile,'file')
%     cfg0 = [];
% else
%     cfg0 = load(cfgFile);
% end

if ~exist('op','var') || isempty(op)
    op = 0;
end

fh = guidata(f);

% new project
if op==0
    preset = fh.preset.Value;
    opts = util.parseParam(preset,0);
    opts.preset = preset;
    
    % read user input
    try
        % if ~strcmp(fh.tmpRes.String,'As preset')
            opts.frameRate = str2double(fh.tmpRes.String);
        % end
        % if ~strcmp(fh.spaRes.String,'As preset')
            opts.spatialRes = str2double(fh.spaRes.String);
        % end
        % if ~strcmp(fh.bdSpa.String,'As preset')
            opts.regMaskGap = str2double(fh.bdSpa.String);
        % end
    catch
        % msgbox('Invalid input');
        % return
    end
    
    try
        pf1 = fh.fIn1.String;
        pf2 = fh.fIn2.String;
        [filepath1,name1,ext1] = fileparts(pf1);
        [filepath2,name2,ext2] = fileparts(pf2);
        f.Name = ['AQUA2: ',name1,' ',name2];
        [datOrg1,datOrg2,opts] = burst.prep1(filepath1,[name1,ext1],filepath2,[name2,ext2],[],opts,ff);
    catch
        msgbox('Fail to load file');
        return
    end
    
    maxPro1 = max(datOrg1,[],3);
    maxPro2 = max(datOrg2,[],3);
    fh.maxPro1 = maxPro1;
    fh.maxPro2 = maxPro2;
    fh.averPro1 = mean(datOrg1,3);
    fh.averPro2 = mean(datOrg2,3);
    fh.showcurves = [];
    opts.maxdF1 = 1;
    opts.maxdF2 = 1;
    guidata(f,fh);
    
    if(isempty(datOrg2))
       opts.singleChannel = true; 
    else
       opts.singleChannel = false; 
    end
    
    opts.alreadyBleachCorrect = 0;
    
    % UI data structure
    [ov,bd,scl,btSt] = ui.proj.prepInitUIStruct(datOrg1,opts); %#ok<ASGLU>
    
    % data and settings
    vBasic = {'opts','scl','btSt','ov','bd','datOrg1','datOrg2'};
    for ii=1:numel(vBasic)
        v0 = vBasic{ii};
        if exist(v0,'var')
            setappdata(f,v0,eval(v0));
        else
            setappdata(f,v0,[]);
        end
    end
    stg = [];
    stg.detect = 0;
end

% read existing project or mat file
if op>0
    if op==1
        fexp = getappdata(f,'fexp');
        tmp = load(fexp);
        res = tmp.res;
        
        % [p00,~,~] = fileparts(fexp);
        % cfg0.outPath = p00;
        % save(cfgFile,'cfg0');
    end
    
    % rescale int8 to [0,1] double
    % dat is for detection, datOrg for viewing
    %res.dat = double(res.dat)/(2^res.opts.bitNum-1);
    if isfield(res,'datOrg1')
        res.datOrg1 = double(res.datOrg1)/(2^res.opts.bitNum-1);
        res.datOrg2 = double(res.datOrg2)/(2^res.opts.bitNum-1);
    else
        res.datOrg1 = double(res.dat1);
        res.datOrg2 = double(res.dat2);
    end
    if isfield(res,'maxVal')
        res.datOrg1 = res.datOrg1*res.maxVal;
        res.datOrg2 = res.datOrg2*res.maxVal;
    end
    
    dat1 = res.datOrg1;
    dat2 = res.datOrg2;
    
    if res.opts.smoXY>0
        for tt=1:size(dat1,3)
            dat1(:,:,tt) = imgaussfilt(dat1(:,:,tt),res.opts.smoXY);
        end
    end
    
    if res.opts.smoXY>0 && ~isempty(dat2)
        for tt=1:size(dat1,3)
            dat1(:,:,tt) = imgaussfilt(dat1(:,:,tt),res.opts.smoXY);
        end
    end
    res.dat1 = dat1;
    res.dat2 = dat2;
    
    waitbar(0.5,ff);
    
    if ~isfield(res,'scl')
        if isfield(res,'bd')
            [~,~,res.scl,res.btSt] = ui.proj.prepInitUIStruct(res.datOrg1,res.opts);
        else
        [~,res.bd,res.scl,res.btSt] = ui.proj.prepInitUIStruct(res.datOrg1,res.opts);
        end
        res.stg = [];
        res.stg.detect = 1;
        res.stg.post = 1;
    else
        [~,~,res.scl,res.btSt] = ui.proj.prepInitUIStruct(res.datOrg1,res.opts,res.btSt);
    end
    
    % reset some settings
    if ~isfield(res,'dbg') || res.dbg==0
        res.btSt.overlayDatSel = 'Events';
    end
    
    opts = res.opts;
    scl = res.scl;
    stg = res.stg;
    ov = res.ov;
    if(~opts.singleChannel)
        fh.nEvt.String = [num2str(numel(res.evt1)),' | ',num2str(numel(res.evt2))];
    else
        fh.nEvt.String = [num2str(numel(res.evt1))];
    end
    
    f.Name = ['AQUA: ',opts.fileName1,' ',opts.fileName2];
    res.btSt.sbs = 0;
    
    fns = fieldnames(res);
    for ii=1:numel(fns)
        f00 = fns{ii};
        setappdata(f,f00,res.(f00));
    end
end

waitbar(1,ff);

% UI
ui.proj.prepInitUI(f,fh,opts,scl,ov,stg,op);

fprintf('Done ...\n');
delete(ff);

end











