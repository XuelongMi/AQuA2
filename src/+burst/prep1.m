function [dat1,dat2,opts] = prep1(p1,f1,p2,f2,rgT,opts,ff)
    %PREP1 load data and estimation noise
    % TODO: use segment by segment processing to reduce the impact of bleaching
    
    bdCrop = opts.regMaskGap;
    
    [filepath1,name1,ext1] = fileparts([p1,filesep,f1]);
    opts.filePath1 = filepath1;
    opts.fileName1 = name1;
    opts.fileType1 = ext1;
    [filepath2,name2,ext2] = fileparts([p2,filesep,f2]);
    opts.filePath2 = filepath2;
    opts.fileName2 = name2;
    opts.fileType2 = ext2;
    
    % read data
    fprintf('Reading data\n');
    if strcmp(ext1,'.mat')
        file = load([p1,filesep,f1]);
        headers = fieldnames(file);
        dat1 = single(file.(headers{1}));
        maxImg1 = -1;
        if(~isempty(f2))
            file = load([p2,filesep,f2]);
            headers = fieldnames(file);
            dat2 = single(file.(headers{1}));
            maxImg2 = -1;
        else
            dat2 = [];
            maxImg2 = [];
        end
    else
        [dat1,maxImg1] = io.readTiffSeq([p1,filesep,f1]);
        if(~isempty(f2))
            [dat2,maxImg2] = io.readTiffSeq([p2,filesep,f2]);
        else
            dat2 = [];
            maxImg2 = [];
        end
    end
    
    if exist('rgT','var') && ~isempty(rgT)
        dat1 = dat1(:,:,rgT);
        if(~isempty(dat2))
            dat2 = dat2(:,:,rgT);
        end
    end
    dat1 = single(dat1);
    dat2 = single(dat2);
    dat1 = dat1(bdCrop+1:end-bdCrop,bdCrop+1:end-bdCrop,:);
    if(~isempty(dat2))
        dat2 = dat2(bdCrop+1:end-bdCrop,bdCrop+1:end-bdCrop,:);
    end
    minDat1 = median(min(dat1,[],[1,2]));
    minDat2 = median(min(dat2,[],[1,2]));
    minDat = min([minDat1,minDat2]);
    dat1 = dat1 - minDat;
    dat2 = dat2 - minDat;
    dat1(dat1<0) = 0;
    dat1(dat2<0) = 0;
    maxDat1 = max(dat1(:)); 
    maxDat2 = max(dat2(:));
    maxDat = max([maxDat1,maxDat2]);
    dat1 = dat1/maxDat;
    dat2 = dat2/maxDat;
    
    if opts.usePG==1
        dat1 = sqrt(dat1);
        dat2 = sqrt(dat2);
    end
    if exist('ff','var')
        waitbar(0.4,ff);
    end
    
    [H,W,T] = size(dat1);
    opts.sz = [H,W,T];
    opts.maxValueDepth = maxImg1;
    opts.maxValueDat = maxDat;
    opts.minValueDat = minDat;
    
end