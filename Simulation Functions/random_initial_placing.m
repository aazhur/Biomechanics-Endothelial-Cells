function [xc,yc] = random_initial_placing(nc,r0,m,r_box)
    xc = zeros(1,nc);
    yc = zeros(1,nc);
    
    R = r_box - r0 - m;
    R2 = R^2;
    L = r_box^2;
    Lc = (m+r0)^2;
    
    while L > R2
        xc(1) = 2*R*rand - R;
        yc(1) = 2*R*rand - R;
        L = xc(1)^2+yc(1)^2; 
    end
    
    for i = 2:nc
        check = 0;
        
        while ~check
            
            xc(i) = 2*R*rand - R;
            yc(i) = 2*R*rand - R;
            L = xc(i)^2+yc(i)^2;
            
            if L < R2 
                
                D = (xc(i)-xc(1:(i-1))).^2+(yc(i)-yc(1:i-1)).^2;    
                
                if min(D) > Lc
                    check = 1;                
                end
                
            end
        end
    end
