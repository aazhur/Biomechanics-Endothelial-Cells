function drawnsave(dir_name,pic_dir_name,itr,effective_r)
    %for i = 2:2
        flnm = [dir_name '\' num2str(itr) '.mat'];
        f = load(flnm, '-mat'); 
        disp(itr);   
                
        figure;    
        hold on;
        axis equal;
        axis off;
        draw_box_flat(f.R,[0,0]);
        hold on;
        draw_cells_alt_flat(f.nc,f.begx(:,1:f.np),f.begy(:,1:f.np),f.endx(:,1:f.np),f.endy(:,1:f.np),f.spx(:,1:f.np),f.spy(:,1:f.np),f.ca(:,1:f.np),f.ba(:,1:f.np),f.xc,f.yc,f.ra,f.rb,f.phi,effective_r);
        
        pic_name = [pic_dir_name '\' num2str(itr) '.png'];
        
        fig = gcf;
        fig.InvertHardcopy = 'off';
        disp(pic_name);
        print(pic_name,'-dpng','-r200');
        
        close(gcf); 
        
    %end   
end