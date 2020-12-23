function [prot_len,endx,endy,endz,spx,spy,spz,ba,es,ca,pa] = prot_attachments(att,nc,np,nd,prot_len,begx,begy,begz,endx,endy,endz,spx,spy,spz,Box_coord,ba,es,ca,pa,ae,xc,yc,zc,ra,rb,phi,R,H,kb,kb_bott,dla)

%% attach/detach
    
    % check cells' detachments
    if sum( sum(ba+ca) ) ~= 0
       lb = (spx-endx).^2+(spy-endy).^2+(spz-endz).^2;
       side = [ones(nc,np),zeros(nc,nd)]; bottom = [zeros(nc,np),ones(nc,nd)]; 
       db = exp(-att.*(lb./(kb.^2.*side+kb_bott.^2.*bottom)));
       p = rand([nc , (np+nd)]) - db;
       p = (p < 0); 
       p(ba == 0 & ca == 0) = 1; %eliminate the possibility of change when not attached
       ba = ba.*p; ca = ca.*p; pa = pa.*p;
       spx = spx.*p + endx.*(1-p); spy = spy.*p + endy.*(1-p); spz = spz.*p + endz.*(1-p);
       es = es.*p + (1-p);
       prot_len = prot_len.*p+dla.*(1-p);
    end  

    %cell clusts
    [Clust,BoxClust] = cluster(xc,yc,zc,prot_len,nc,R,ra,rb); 
    
    % Now check for attachments
    Endx_d = endx; Endy_d = endy; Endz_d = endz; 
    % make impossible for not growing detached portrusions to be found as
    % attached, as well as recheck already attached once
    Endx_d(ca ~= 0 | ba ~= 0 | es == 0 ) = 3*R; 
    Endy_d(ca ~= 0 | ba ~= 0 | es == 0 ) = 3*R;
    Endz_d(ca ~= 0 | ba ~= 0 | es == 0 ) = 3*R;

    for i=1:nc
   
        if sum(Clust(i,:)) ~= 0
            
            %disp(size(Clust(i,:)));
            idx = find(Clust(i,:));
            Xc = xc(idx); Yc = yc(idx); Zc = zc(idx);
            Ra = ra(idx); Rb = rb(idx); Phi = phi(idx);
                
            % check every protrusion if inside another cell
            [start_in,start_on] = inellipse(np,begx(i,1:np),begy(i,1:np),numel(idx),Xc,Yc,Ra,Rb,Phi);
            start = start_in+start_on;
            [in,on] = inellipse(np,Endx_d(i,1:np),Endy_d(i,1:np),numel(idx),Xc,Yc,Ra,Rb,Phi);
            in = in.*(start == 0); on = on.*(start == 0);
            [~,n] = size(idx);
            idx = diag(idx)*ones(n,np);
            where = ones(n,np)*diag(ae(i,:));
            idx = idx.*(idx == where);
            
            if n == 1
                att = (in+on)'.*idx;
            else
                att = (in+on)'.*idx;
                att = max(att);
            end
            
            if sum(att) ~= 0
                
                num_att = find(att > 0);
                att_ = att(att > 0);
                n = zeros(1,(np+nd));

                for j = 1:numel(att_)
                    %if ismember(i,ca(att_(j),:)) == 0
                       L_end = (Endx_d(i,num_att(j))-begx(att_(j),1:np)).^2+(Endy_d(i,num_att(j))-begy(att_(j),1:np)).^2;
                       L_beg = (begx(i,num_att(j))-begx(att_(j),1:np)).^2+(begy(i,num_att(j))-begy(att_(j),1:np)).^2;
                       [~,n(num_att(j))] = min(L_end+L_beg); 
                    %else
                        %n(num_att(j)) = -1;
                    %end   
                end
                pa(i,1:np) = pa(i,1:np) + n(1:np).*(n(1:np) > 0);
                
                spx_add = zeros(1,(np+nd)); spy_add = zeros(1,(np+nd)); spz_add = zeros(1,(np+nd));
                
                for j = 1:numel(att_)
                    %if n(num_att(j)) > 0
                        nca = att(num_att(j)); npa = n(num_att(j));
                        spx_add(num_att(j)) = begx(nca,npa);
                        spy_add(num_att(j)) = begy(nca,npa); 
                        spz_add(num_att(j)) = begz(nca,npa);
                    %end    
                end    
                ca(i,1:np) = ca(i,1:np) + att.*(ca(i,1:np) == 0).*(n(1:np) > 0);
                es(i,1:np) = es(i,1:np).*(att == 0); 
                
                spx(i,:) = spx(i,:).*(n == 0)+begx(i,:).*(n == -1)+spx_add; 
                spy(i,:) = spy(i,:).*(n == 0)+begy(i,:).*(n == -1)+spy_add; 
                spz(i,:) = spz(i,:).*(n == 0)+begz(i,:).*(n == -1)+spz_add;
                endx(i,:) = endx(i,:).*(n == 0)+begx(i,:).*(n == -1)+spx_add; 
                endy(i,:) = endy(i,:).*(n == 0)+begy(i,:).*(n == -1)+spy_add;
                endz(i,:) = endz(i,:).*(n == 0)+begz(i,:).*(n == -1)+spz_add;
                prot_len(i,:) = sqrt((begx(i,:)-endx(i,:)).^2+(begy(i,:)-endy(i,:)).^2+(begz(i,:)-endz(i,:)).^2);
            end
        end 
    end 
    

    Endx_d = endx; Endy_d = endy; Endz_d = endz;
    % zero, so cannot be wrongly found attached to a box, when in search or
    % attached to a cell
    Endx_d(ca ~= 0 | ba ~= 0 | es == 0) = 0;
    Endy_d(ca ~= 0 | ba ~= 0 | es == 0) = 0;
    Endz_d(ca ~= 0 | ba ~= 0 | es == 0) = R;
           
    % check every protrusion if stuck to the box outer rim
    for i = 1:nc
        if BoxClust(i) == 1            
            
            [in,on] = inellipse(np,Endx_d(i,1:np),Endy_d(i,1:np),1,0,0,R,R,0);
            att = (~in')+on';
        
            if sum(att) ~= 0
                for j = 1:numel(att)
                    if att(j) == 1
                          L = sqrt((endx(i,j)-Box_coord(:,1)).^2+(endy(i,j)-Box_coord(:,2)).^2); 
                          M = max(L);
                          L(Box_coord(:,3) ~= H) = M;
                          [~,n] = min(L); 
                          spx(i,j) = Box_coord(n,1); spy(i,j) = Box_coord(n,2); spz(i,j) = H;
                          endx(i,j) = Box_coord(n,1); endy(i,j) = Box_coord(n,2); endz(i,j) = H;
                          prot_len(i,j) = sqrt((begx(i,j)-endx(i,j))^2+(begy(i,j)-endy(i,j))^2+(begz(i,j)-endz(i,j))^2);
                          ba(i,j) = n; es(i,j) = 0;
                    end    
                end    
            end
        end
    end
 
    %check every protrusion stuck to the bottom of the box
    for i = 1:nc           
        att = Endz_d(i,(np+1):end) <= 0;

        if sum(att) ~= 0
            for j = 1:numel(att)
                if att(j) == 1
                      k = j + np;
                      L = sqrt((endx(i,k)-Box_coord(:,1)).^2+(endy(i,k)-Box_coord(:,2)).^2+(endz(i,k)-Box_coord(:,3)).^2); 
                      [~,n] = min(L); 
                      spx(i,k) = Box_coord(n,1); spy(i,k) = Box_coord(n,2); spz(i,k) = Box_coord(n,3);
                      endx(i,k) = Box_coord(n,1); endy(i,k) = Box_coord(n,2); endz(i,k) = Box_coord(n,3);
                      prot_len(i,k) = sqrt((begx(i,k)-endx(i,k))^2+(begy(i,k)-endy(i,k))^2+(begz(i,k)-endz(i,k))^2);
                      ba(i,k) = n; es(i,k) = 1;
                end    
            end    
        end
    end
  
end   
    
% check what cells and whether the box are in vicinity of each cell
function [Clust, BoxClust] = cluster(xc,yc,zc,prot_len,nc,R,ra,rb)
    r = max(prot_len,[],2)'+0.5*(ra+rb+abs(ra-rb));

    radius = sqrt(xc.^2 + yc.^2)+r;
    Dist = R - radius;
    BoxClust = (Dist <= 0);

    r = diag(r)*ones(nc,nc)+ones(nc,nc)*diag(r);
    X = (diag(xc)*ones(nc,nc)-ones(nc,nc)*diag(xc)).^2;
    Y = (diag(yc)*ones(nc,nc)-ones(nc,nc)*diag(yc)).^2; 
    Z = (diag(zc)*ones(nc,nc)-ones(nc,nc)*diag(zc)).^2; 
    Intersec = sqrt(X+Y+Z) - r;
    Clust = (Intersec <= 0) ;
    Clust(logical(eye(size(Clust)))) = 0;
end