function [pdf,binSize] = getDistribution(difdF,binSize)
    %% here assume the distribution is symmetric
    difdF = abs(difdF(:));
    if(~exist('binSize','var')||isempty(binSize))
        SampleSizeRatio = (length(difdF))^(1/3);
        binSize = 1.5*2*max(difdF)/SampleSizeRatio;
    end
    
    N0 = ceil((max(difdF)-0.5*binSize)/binSize);    % number of half side
    edges = [0,([1:N0+1]-0.5)*binSize]; % set edges
    pdf_half = histcounts(difdF,edges);
    pdf_half(1) = pdf_half(1)*2 - sum(difdF==0);  % avoid count 0 twice
    pdf = [fliplr(pdf_half(2:end)),pdf_half];
    pdf = pdf/sum(pdf);
    pdf = pdf';
    if(mod(N0,2)~=0)
        pdf = [0;pdf;0];
    end
end