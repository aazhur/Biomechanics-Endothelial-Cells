function M = box_surface(R,H,center)
    Rx = R; Ry = R;
    Hz = H;

    M = zeros(360*200,3);
    theta = linspace(0,360,360);
    height = linspace(0,Hz,100);
    radius = linspace(0,Rx,100);

    for i = 0:99
        M(360*i+1:360*(i+1),1) = Rx*cosd(theta)+center(1);
        M(360*i+1:360*(i+1),2) = Ry*sind(theta)+center(2);
        M(360*i+1:360*(i+1),3) = height(i+1);
    end

    for i = 100:199
        M(360*i+1:360*(i+1),1) = radius(i-99)*cosd(theta)+center(1);
        M(360*i+1:360*(i+1),2) = radius(i-99)*sind(theta)+center(2);
    end
end