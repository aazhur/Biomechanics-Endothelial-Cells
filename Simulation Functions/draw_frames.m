clear;
prompt = 'Enter directory name you would like to view: ';
dir_name = uigetdir(cd,prompt);

test_info = [dir_name '\test_info.mat']; %Windows
load(test_info, '-mat');

frm = Nitr/ditr + 1;
gx = 10.0*ones(1,nc);
disp(gx);
run = 1;

prompt = 'Enter directory name to save pictures: ';
pic_dir_name = input(prompt,'s');
s = mkdir('Y:\tsygankov-lab\Nastya Zhurikhina\Simulations\', pic_dir_name);
%%Windows

for itr = frm:frm
    
        disp(itr);
        drawnsave(dir_name,pic_dir_name,itr,run,gx)
end  