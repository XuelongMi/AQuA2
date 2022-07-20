function [dF,xBias] = getDfBlk_linear_removeNaN2(datIn,cut,movAvgWin,nonvalid)
    
    % non-valid is to show the region caused by registration edges
    % grow 5 time points for nan valid part in temporal

    [H,W,T] = size(datIn);
    datMA = datIn;
    datMA(nonvalid) = nan;
    datMA = movmean(datMA,movAvgWin,3,'omitnan');
    datMA(nonvalid) = nan;
    step = round(0.5*cut);
    
    maxV = max(datIn(:));
    dF = zeros(size(datIn),'single');
    nSegment = ceil(T/step)-1;
    for k = 1:nSegment
         t0 = 1 + (k-1)*step;
         t1 = min(T,t0+cut);
         [curV,curP] = nanmin(datMA(:,:,t0:t1),[],3);
         curV(isnan(curV)) = 0;
         se = offsetstrel('ball',5,5);
         maxCurV = imdilate(curV,se) - 5;
         maxCurV = imgaussfilt(maxCurV,2);
%         maxCurV = curV;%imgaussfilt(curV,2);
%          maxCurV = zeros(size(curV));
%          for x = 1:H
%              for y = 1:W
%                  xRange = max(x-dist,1):min(x+dist,H);
%                  yRange = max(y-dist,1):min(y+dist,W);
%                  maxCurV(x,y) = max(max(curV(xRange,yRange)));
%              end 
%          end
         curP = curP + t0 - 1;

         for x = 1:H
             for y = 1:W
                 if(k==1)
                     dF(x,y,1:curP(x,y)) = maxCurV(x,y);
                 else
                     mt1 = preP(x,y);
                     mt2 = curP(x,y);
                    if(mt1<mt2)
                        value1 = preV(x,y);
                        value2 = maxCurV(x,y);
                        dF(x,y,mt1:mt2) = value1 + (value2-value1)/(mt2-mt1)*[0:mt2-mt1]; 
                    end
                 end
                 
                 if(k==nSegment)
                     dF(x,y,curP(x,y):T) = maxCurV(x,y);
                 end
             end
         end
         preP = curP;
         preV = maxCurV;
    end
    
    xx = randn(10000,cut); %% bias need to be reconsidered
    xxMA = movmean(xx,movAvgWin,2);
    xxMin = min(xxMA,[],2);
    xBias = nanmean(xxMin(:));
    clear datMA;
    dF = datIn - dF;    % for save memory, dF is F0 before this line
    dF(nonvalid) = 0;
end