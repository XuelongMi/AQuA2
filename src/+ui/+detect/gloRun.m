function gloRun(~,~,f)
% z scores and filtering

fh = guidata(f);
opts.detectGlo = fh.detectGlo.Value;
if(opts.detectGlo)
    fprintf('Detecting global signal ...\n')
    ff = waitbar(0,'Merging ...');
    delete(ff);
    fprintf('Done\n')
end
end