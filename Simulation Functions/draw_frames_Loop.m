%function draw_frames_Loop    
    close all;
    last = 3;
    %prompt = 'Enter directory name you would like to view: ';
    %folder_name = uigetdir(cd,prompt)
    for n2 = 20
        n1 = 6000; n5 = 100; n3 = 0.5; n4 = 120; n6 = 3;
        dir_name = strcat('Y:\tsygankov-lab\Nastya Zhurikhina\Simulations\Testing 3d new scaling\08.10.18\wt_kb=', num2str(n1), ';kb_bott=', num2str(n2), ';ks=', num2str(n3), ';r=', num2str(n4), ';kz=', num2str(n5), ';ka=', num2str(n6));
        test_info = [dir_name '\test_info.mat']; %Windows
        load(test_info, '-mat');

        frm = Nitr/ditr + 1;

        pic_dir_name = strcat('Y:\tsygankov-lab\Nastya Zhurikhina\Simulations\Pictures for Testing 3d new scaling\08.10.18\wt_r=', num2str(n1), ';kb_bott=', num2str(n2), ';ks=', num2str(n3), ';r=', num2str(n4), ';kz=', num2str(n5), ';ka=', num2str(n6));
        mkdir(pic_dir_name);

        parfor itr = 1:frm
                %try
                    drawnsave(dir_name,pic_dir_name,itr,120);  
                %catch
                %end
        end 
        
    end 