function draw_cells_alt_flat(nc,begx,begy,endx,endy,spx,spy,ca,ba,xc,yc,ra,rb,phi,effective_r)
    
    is_att = ca+ba; is_att(is_att ~= 0) = 1;

    begx_a = (begx(is_att == 1))'; begy_a = (begy(is_att == 1))';
    spx = (spx(is_att == 1))'; spy = (spy(is_att == 1))';

    x_prot = [begx_a; spx];
    y_prot = [begy_a; spy];

    x_notatt = [(begx(is_att == 0))'; (endx(is_att == 0))'];
    y_notatt = [(begy(is_att == 0))'; (endy(is_att == 0))'];
    
    %att_color = [0.55,0.0706,0.1217];
    protrusion = plot(x_prot, y_prot, 'k', 'LineWidth', 1);
    %set(protrusion, 'Color', att_color);

    %{
    prot_color = [0,0.5,0];
    protrusion = plot(x_notatt, y_notatt, 'k', 'LineWidth', 0.5);
    %set(protrusion, 'Color', prot_color);
   %}
    [Xo, Yo] = ellipse(xc, yc, ra, rb, phi, nc);
    for i=1:nc
            fill(Xo(i,:), Yo(i,:), 'k');
            %plot(xc(i) + effective_r*cosd(0:10:360), yc(i) + effective_r*sind(0:10:360), '-b');  
    end
  
end

function [Xo, Yo] = ellipse(x0, y0, ra, rb, phi, nc)
         teta = 0:360;
 
         xo = ra'*cosd(teta); yo = rb'*sind(teta);
         Xo = diag(x0)*ones(nc , 361) + diag(cosd(phi))*xo - diag(sind(phi))*yo;
         Yo = diag(y0)*ones(nc , 361) + diag(sind(phi))*xo + diag(cosd(phi))*yo;
end