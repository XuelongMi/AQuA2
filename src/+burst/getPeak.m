function tw = getPeak(curve,rgT,s1)

    peak0 = max(curve(rgT));
    it0 = min(rgT);
    it1 = max(rgT);
    thr = 3*s1;
    cnt = 1;
    tw = cell(0);
    while 1
       [peak,its] = nanmax(curve(rgT)); 
       its = its + min(rgT)-1;
        if(isempty(peak)||isnan(peak)||peak<0.8*peak0)
            break;
        end
        
        %% seed grow, left
        t0 = it0;
        base0 = peak;
        st0 = 0;
        for tt = its-1:-1:it0
           if curve(tt)<base0
                base0 = curve(tt);
                iMin = tt;
            end
            if curve(tt)>peak0
                t0 = iMin;
                break
            end
            if peak0-curve(tt)>thr
                st0 = 1;
            end
            if curve(tt)==-100 || (curve(tt)-base0>thr && st0==1)
                t0 = iMin;
                break
            end
            if isnan(curve(tt))
                t0 = tt;
                break
            end 
        end
        
        %% seed grow, right
        st0 = 0;
        t1 = it1;
        base1 = peak;
        for tt = its+1:it1
           if curve(tt)<base1
                base1 = curve(tt);
                iMin = tt;
            end
            if curve(tt)>peak0
                t1 = iMin;
                break
            end
            if peak0-curve(tt)>thr
                st0 = 1;
            end
            if curve(tt)==-100 || (curve(tt)-base1>thr && st0==1)
                t1 = iMin;
                break
            end
            if isnan(curve(tt))
                t1 = tt;
                break
            end 
        end
        
        rgT0 = t0:t1;
        tw{cnt} = rgT0;
        %% update
        curve(rgT0) = -100;
        cnt = cnt+1;
    end


end