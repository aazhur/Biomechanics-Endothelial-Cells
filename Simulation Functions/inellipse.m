function [in,on] = inellipse(sp,xp,yp,cp,xc,yc,ra,rb,phi)
    
    rb = (diag(rb)*ones(cp,sp))';
    ra = (diag(ra)*ones(cp,sp))';
    Cos = (diag(cosd(phi))*ones(cp,sp))';
    Sin = (diag(sind(phi))*ones(cp,sp))';
    X = diag(xp)*ones(sp,cp)-(diag(xc)*ones(cp,sp))';
    Y = diag(yp)*ones(sp,cp)-(diag(yc)*ones(cp,sp))';
    
    L = (Cos.*X+Sin.*Y).^2./ra.^2 + (Sin.*X-Cos.*Y).^2./rb.^2;
    in = L < 1;
    on = L == 1;
        
end