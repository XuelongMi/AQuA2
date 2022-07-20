function df0ip = imputeMov_Fast(df0,validMap)

[H,W,T] = size(df0);
df0ip = reshape(df0,[],T);
pix = find(validMap);
for ii=1:numel(pix)
    curPix = pix(ii);
    x0 = df0ip(curPix,:);
    for tt=2:T
        if isnan(x0(tt))
            x0(tt) = x0(tt-1);
        end
    end
    for tt=T-1:-1:1
        if isnan(x0(tt))
            x0(tt) = x0(tt+1);
        end
    end
    df0ip(curPix,:) = x0;
end
df0ip = reshape(df0ip,[H,W,T]);
df0ip(isnan(df0ip)) = 0;

end