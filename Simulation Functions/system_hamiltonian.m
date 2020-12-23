function H = system_hamiltonian(v,k,nc,np,nd,prot_len,begx,begy,begz,endx,endy,endz,spx,spy,spz,ba,pa,ca,xc,yc,zc,ra,rb,phi,shift,V,N,ka,kz,ks,q)
    %define which coordinate to change
    if v == 1
        xc(k) = xc(k) + q;
    elseif v == 2
        yc(k) = yc(k) + q;
    elseif v == 3
        ra(k) = ra(k) + q;
    elseif v == 4
        rb(k) = rb(k) + q;
    elseif v == 5
        phi(k) = phi(k) + q;
    elseif v == 6
        shift(k) = shift(k) + q;    
    end
    
    % contribution from cell deformation
    rc = (3*V./(ra.*rb*4*pi));
    R = (3*V/(4*pi)).^(1/3);
    E = 0.5*(ka*(ra-R).^2+ka*(rb-R).^2+kz*(rc-R).^2);
    H0 = sum(E);
    %{
    if(v == 3)
        disp('H0');
        disp([H0,q]);
    end    
    %}
    % update coordinates
    [endx,endy,endz,spx,spy,spz] = one_cell_change(v,k,q,np,nd,prot_len,begx,begy,begz,endx,endy,endz,spx,spy,spz,pa,ca,ba,xc(k),yc(k),zc(k),ra(k),rb(k),phi(k),shift(k),N);
    % index of attached protrusions
    side = [ones(nc,np),zeros(nc,nd)]; bottom = [zeros(nc,np),ones(nc,nd)]; 	
    E = 0.5*ks*((endx-spx).^2+(endy-spy).^2+(endz-spz).^2).*side+0.5*0.5*((endx-spx).^2+(endy-spy).^2+(endz-spz).^2).*bottom;
    H = H0+sum(sum(E));
    %{
    if(v == 3)
        disp('H');
        disp([H,q]);
    end 
    %}
    
