function getOutputFolder(~,~,f)
    
    opts = getappdata(f,'opts');
    
    % SP, 18.07.16
    definput = {'_AQuA2'};
    selname = inputdlg('Type desired suffix for output files name:',...
        'Output files',[1 75],definput);
    
    selname = char(selname);
    if isempty(selname)
        selname = '_AQuA2';
    end
    file0 = [opts.fileName1,selname];
    clear definput selname
    
    %file0 = [opts.fileName,'_AQuA']; SP, 18.07.16
    selpath = uigetdir(opts.filePath1,'Choose output folder');
    path0 = [selpath,filesep,opts.fileName1];
    if ~isnumeric(selpath)
        ui.proj.saveExp([],[],f,file0,path0);
    end
    
end