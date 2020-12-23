function draw_box_flat(R,center)

    teta = linspace(0,360,360);

    xo = R*cosd(teta)+center(1); 
    yo = R*sind(teta)+center(2);

    plot(xo,yo,'w');
    fill(xo, yo, 'w','LineStyle','none');
    axis equal;
         
         %{
         gb = 1 - sqrt(sum(ba(:))/(np*nc));
         edge_color = [1,0.89*gb,0.8*gb];
         edge_c = plot(x_box, y_box, 'Linewidth', 3.0);
         set(edge_c,'Color',edge_color);
         %}
end