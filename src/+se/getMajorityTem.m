function TW = getMajorityTem(curve,TW,pix,sz,forMerging)

    if(~exist('forMerging','var'))
        forMerging = 0;
    end

    t00 = min(TW);
    t11 = max(TW);
    ext = 3;
    ext00 = 0;
    H = sz(1); W = sz(2); T = sz(3);
    [~,~,it]  = ind2sub([H,W,T],pix);
    t0 = max(1,min(it)-ext00); t1 = min(T,max(it)+ext00);

    %% Spatial majority
    s0 = sqrt(median((curve(2:end)-curve(1:end-1)).^2)/0.9133);
    [~,tPeak] = max(curve(t00:t11));
    tPeak = tPeak + t00 - 1;

    % start time
    if(forMerging)
        [minV,tw0] = min(curve(t00:max(t00,tPeak-ext)));
        tw0 = tw0 + t00 - 1;
    else
        % new change here
        tw0 = t00; minV = curve(tw0);
    end
    ts = tw0;
    
    for t = tw0:-1:t0
        if(curve(t)<minV)
            minV = curve(t);
            ts = t;
        else
            if(curve(t)-minV>=3*s0)
                break;
            end
        end
    end

    % end time
    if(forMerging)
        tRs = min(t11,tPeak+ext);
        [minV,tw1] = min(curve(tRs:t11));
        tw1 = tw1 + tRs - 1;
    else
        % new change here
        tw1 = t11;minV = curve(tw1);
    end
    te = tw1;
    for t = tw1:t1
        if(curve(t)<minV)
            minV = curve(t);
            te = t;
        else
            if(curve(t)-minV>=3*s0)
                break;
            end
        end
    end
       
%     %% let it symmetric
%     if(curve(ts)>curve(te))
%         while(curve(te)<curve(ts))
%             te = te-1;
%         end
%     else
%         while(curve(te)>curve(ts))
%             ts = ts + 1;
%         end
%     end
%     ts = min(ts,TW(1));
%     te = max(te,TW(end));
    TW = ts:te;
end