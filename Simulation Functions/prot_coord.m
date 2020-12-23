function [begx, begy, begz, endx, endy, endz, spx, spy, spz, L] = prot_coord(np,nc,nd,N,prot_len,Box_coord,spx,spy,spz,ba,ca,pa,xc,yc,zc,ra,rb,phi,shift)
        
        L = zeros(nc,(np+nd));
        begx = zeros(nc,(np+nd)); begy = zeros(nc,(np+nd)); begz = zeros(nc,(np+nd));
        endx = zeros(nc,(np+nd)); endy = zeros(nc,(np+nd)); endz = zeros(nc,(np+nd));

        %define protusions beg and end
        alpha = 360*(1:np)/np;

        % begining
        Xbeg = diag(ra)*cosd(diag(shift)*ones(nc,np)+(ones(nc,1))*alpha);
        Ybeg = diag(rb)*sind(diag(shift)*ones(nc,np)+(ones(nc,1))*alpha);
      
        %ending when no attached
        Xend = Xbeg + prot_len(:,1:np).*cosd(diag(shift)*ones(nc,np)+(ones(nc,1))*alpha);
        Yend = Ybeg + prot_len(:,1:np).*sind(diag(shift)*ones(nc,np)+(ones(nc,1))*alpha);
            
        %begining
        begx(:,1:np) = diag(xc)*ones(nc,np) + diag(cosd(phi))*Xbeg - diag(sind(phi))*Ybeg;
        begy(:,1:np) = diag(yc)*ones(nc,np) + diag(sind(phi))*Xbeg + diag(cosd(phi))*Ybeg;
        begz(:,1:np) = diag(zc)*ones(nc,np);
        
        for i = 1:nc
            [x,y,~] = StartingPoints(N,0.9*ra(i),0.9*rb(i));
            begx(i,(np+1):end) = xc(i) + x*cosd(phi(i)) - y*sind(phi(i)); 
            begy(i,(np+1):end) = yc(i) + x*sind(phi(i)) + y*cosd(phi(i)); 
            begz(i,(np+1):end) = zc(i);
        end  
        
        %ending
        endx(:,1:np) = diag(xc)*ones(nc,np) + diag(cosd(phi))*Xend - diag(sind(phi))*Yend;
        endy(:,1:np) = diag(yc)*ones(nc,np) + diag(sind(phi))*Xend + diag(cosd(phi))*Yend;
        endz(:,1:np) = diag(zc)*ones(nc,np);  
    
        endx(:,(np+1):end) = begx(:,(np+1):end);
        endy(:,(np+1):end) = begy(:,(np+1):end);
        endz(:,(np+1):end) = begz(:,(np+1):end) - prot_len(:,(np+1):end);
           
    parfor i=1:nc
        for j=1:(np+nd)
            if ca(i,j) ~= 0 || ba(i,j) ~= 0                
                %check for changes in "sticking" coordinates
                if ca(i,j) ~= 0
                    k = ca(i,j); n = pa(i,j);
                    spx(i,j) = begx(k,n); spy(i,j) = begy(k,n); spz(i,j) = begz(k,n);
                elseif ba(i,j) ~= 0 && j <= np
                    n = ba(i,j);
                    spx(i,j) = Box_coord(n,1); spy(i,j) = Box_coord(n,2); spz(i,j) = Box_coord(n,3);    
                end    
                L(i,j) = sqrt((begx(i,j) - spx(i,j))^2 + (begy(i,j) - spy(i,j))^2 + (begz(i,j) - spz(i,j))^2);

                endx(i,j) = begx(i,j) + prot_len(i,j)*(spx(i,j) - begx(i,j))/L(i,j);
                endy(i,j) = begy(i,j) + prot_len(i,j)*(spy(i,j) - begy(i,j))/L(i,j); 
                endz(i,j) = begz(i,j) + prot_len(i,j)*(spz(i,j) - begz(i,j))/L(i,j); 
            else
                spx(i,j) = endx(i,j);
                spy(i,j) = endy(i,j);
                spz(i,j) = endz(i,j);
            end

        end
    end 
    
end
