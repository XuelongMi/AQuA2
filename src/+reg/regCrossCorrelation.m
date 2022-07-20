function [data1,data2,tforms] = regCrossCorrelation(data1,data2)
    ref = mean(data1(:,:,1:10),3);
    ref = median(ref(:)) - ref; % align dark
    [H,W,T] = size(data1);
    tforms = cell(T,1);
    x0 = 0;
    x1 = 0;
    y0 = 0;
    y1 = 0;
    
    % concise cross correlation
    for t = 2:T
        moving = data1(:,:,t);
        moving = median(moving(:)) - moving;
        matrix = calCC(moving,ref);
        [~,id] = max(matrix(:));
        [hShift,wShift] = ind2sub(size(matrix),id); 
        tform = eye(3);
        tform(3,1) = W-wShift;
        tform(3,2) = H-hShift;
        x0 = min(W-wShift,x0);
        x1 = max(W-wShift,x1);
        y0 = min(H-hShift,y0);
        y1 = max(H-hShift,y1);
        tforms{t} = affine2d(tform);
    end

    for t = 2:T 
        data1(:,:,t) = imwarp(data1(:,:,t),tforms{t},'OutputView',imref2d([H,W]));
    end
    data1 = data1(y1+1:end+y0,x1+1:end+x0,:);
    
    if(~isempty(data2))
        for t = 2:T 
            data2(:,:,t) = imwarp(data2(:,:,t),tforms{t},'OutputView',imref2d([H,W]));
        end
        data2 = data2(y1+1:end+y0,x1+1:end+x0,:);
    end
end