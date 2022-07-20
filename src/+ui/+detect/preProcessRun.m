function preProcessRun(~,~,f)
% active voxels detection and update overlay map

opts = getappdata(f,'opts');

preSetting = getappdata(f,'preSetting');
fh = guidata(f);
if(isempty(preSetting) || ~isfield(opts,'alreadyProprecess') || ~opts.alreadyProprecess || ...
        preSetting.registrateCorrect~=fh.registrateCorrect.Value || preSetting.bleachCorrect~=fh.bleachCorrect.Value)
    % 'if' is to judge whether this step is already done. Since this step
    % is time-consuming.
    opts = getappdata(f,'opts');
    datOrg1 = getappdata(f,'datOrg1');
    datOrg2 = getappdata(f,'datOrg2');
    ff = waitbar(0,'Correcting ...');
    preSetting.registrateCorrect = fh.registrateCorrect.Value;
    preSetting.bleachCorrect = fh.bleachCorrect.Value;
    setappdata(f,'preSetting',preSetting);
    if(fh.registrateCorrect.Value == 2)
        [datOrg1,datOrg2] = reg.regCrossCorrelation(datOrg1,datOrg2);
    elseif(fh.registrateCorrect.Value == 3)
        if(opts.singleChannel)
            [datOrg1,datOrg2] = reg.regCrossCorrelation(datOrg1,datOrg2);
        else
            [datOrg2,datOrg1] = reg.regCrossCorrelation(datOrg2,datOrg1);
        end
    end

    if(fh.bleachCorrect.Value==2)
        [datOrg1] = burst.bleach_correct(datOrg1);
        if(~opts.singleChannel)
            [datOrg2] = burst.bleach_correct(datOrg2);
        end
    elseif(fh.bleachCorrect.Value==3)
        [datOrg1] = burst.bleach_correct2(datOrg1,opts);
        if(~opts.singleChannel)
            [datOrg2] = burst.bleach_correct2(datOrg2,opts);
        end
    end

    opts.bleachCorrect = fh.bleachCorrect.Value;
    opts.registrateCorrect = fh.registrateCorrect.Value;
    opts.alreadyProprecess = true;
    opts.sz = size(datOrg1);
    fh.averPro1 = mean(datOrg1,3);
    fh.averPro2 = mean(datOrg2,3);
    guidata(f,fh);
    setappdata(f,'datCorrect1',datOrg1);
    setappdata(f,'datCorrect2',datOrg2);
    setappdata(f,'opts',opts);
    waitbar(1,ff);
    delete(ff);
    fprintf('Done\n');
end

end