function [x_dir,y_dir,z_dir] = dirGenerate(conn)
    if(conn==6)
        x_dir = [-1;1;0;0;0;0];
        y_dir = [0;0;-1;1;0;0];
        z_dir = [0;0;0;0;-1;1];
    else
        x_dir = zeros(27,1);
        y_dir = zeros(27,1);
        z_dir = zeros(27,1);
        cnt = 1;
        
        for x = -1:1
            for y = -1:1
                for z = -1:1
                    x_dir(cnt) = x;
                    y_dir(cnt) = y;
                    z_dir(cnt) = z;
                    cnt = cnt + 1;
                end
            end
        end
        
        x_dir(14) = [];
        y_dir(14) = [];
        z_dir(14) = [];
    end

end