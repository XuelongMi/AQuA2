function [dFResize,H_resize,W_resize] = normalizeAndResize(dFOrg,opts)
    scaleRatios = opts.scaleRatios;
    dFResize = cell(numel(scaleRatios),1);
    [H,W,T] = size(dFOrg);
    H_resize = zeros(numel(scaleRatios),1);
    W_resize = zeros(numel(scaleRatios),1);
     
    %% Rescl dFOrg
    for i = 1:numel(scaleRatios)
        scaleRatio = scaleRatios(i);
        dFONN = se.myResize(dFOrg,1/scaleRatio);
        H0 = floor(H/scaleRatio);
        W0 = floor(W/scaleRatio);
        dFONN = dFONN(1:H0,1:W0,:);
        H_resize(i) = H0;
        W_resize(i) = W0;
        dFResize{i} = dFONN*scaleRatio;
    end
end