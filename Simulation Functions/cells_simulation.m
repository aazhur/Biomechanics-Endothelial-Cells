function cells_simulation(kb,kb_bott,kz,ka,ks,effective_r,close_r,effective_n,close_n,type)
%% variables and parameters used
% nc - number of cells
% xc - cell centers' x_coord
% yc - cell centers' y_coord
% ra - main elliptic diameter
% S0 - initial mean area2
% S - area of cell
% phi - angel of cell rotation
% r0 - initial ellipse radius
% Xo - matrix of cells' exterior x-coordinates
% Yo - matrix of cells' exterior x-coordinates

% np - number of protrusions
% prot_len - lengths and attachment status of every protrusion
% ca - cell attachments of protrusions
% ae - aspirational elongation
% pa - number of cell's protrusion to wich a corresponding protsrusion is
% attached
% spx - end x coordinates of springs
% spy - end y coordinates of springs
% dl - portrusions elongation/reduction value
% l0 - initial portrusions length
% k - random generator coefficient

% q - dx for Hamilton
% ka - a constant factor characteristic of ellipse
% ks - a constant factor characteristic of the spring

%make each simulation unique
rng shuffle;

%cell-ecm constant values
%% box
edge_dist = 50;
coord_len = 863 + edge_dist;
%coord_len = 154 + edge_dist;

%% set box
teta = 0:360;
R = coord_len - edge_dist;
%effective_r = 55;
disp(R);
H = 2.6; %2*mld vlaue (mean protrusions length when detached)
Box_coord = box_surface(R,H,[0,0]);

%% main parameters of ellipse
np = ceil(pi*effective_r/sqrt(40));
nc = 466;
%nc = 15;
N = 4;
V0 = ((sqrt(40))^3)*4*pi/3;
V = V0*ones(1 , nc);
ra = sqrt(40)*ones(1 , nc);
rb = sqrt(40)*ones(1 , nc);
rc = sqrt(40)*ones(1 , nc);
shift = zeros(1, nc);
r0 = sqrt(40);
[~,~,nd] = StartingPoints(N,0.9*r0,0.9*r0);
dl = 1;
m = 2*r0; %minimal initial cell distance
[xc,yc] = random_initial_placing(nc,r0,m,R);
zc = H*ones(1,nc);

%% set cells turned to the closest neighbors
phi = angle_to_neighbor(xc,yc,nc,R);

%% Hamilton gradient
q = 0.01;
%ka = 5.0;
%kz = 600.0;
%ks = 0.5;
t = 0.025;

%{
%% protrusions length parameters
ua = 0.01;   % (mean-length) rate of switching to "growth" mode (when attached)
wa = 0.99;   % (mean-length) rate of switching to "shrink" mode (when attached)
    
ud = 0.5;   % (mean-length) rate of switching to "growth" mode (when dettached)
wd = 0.5;   % (mean-length) rate of switching to "shrink" mode (when dettached)
    
mld = ml*[ones(nc,np),zeros(nc,nd)] + 2.6*[zeros(nc,np),ones(nc,nd)];  %mean length when detached; 
mla = 0;  %mean length when attached
%}
%% portrusions main parameters
k = 0.1;
kl = 3.5;
l0 = 0*[ones(nc,np),zeros(nc,nd)] + 3*[zeros(nc,np),ones(nc,nd)];
dla = 0.21; % dl when attached
% set initial protrusions length and attachments
prot_len = l0.*ones(nc , (np + nd));

%% springs
%kb = 3300;
%kb_bott = kb_bott*ones(1,nc);
att = (type > 0)*[ones(nc,np),zeros(nc,nd)] + 1*[zeros(nc,np),ones(nc,nd)];
spx = zeros(nc , (np+nd));
spy = zeros(nc , (np+nd));
spz = zeros(nc , (np+nd));
ca = zeros(nc , (np+nd));
ba = zeros(nc , (np+nd));
pa = zeros(nc , (np+nd));
es = ones(nc , (np+nd));

[begx, begy, begz, endx, endy, endz, spx, spy, spz, L] = prot_coord(np,nc,nd,N,prot_len,Box_coord,spx,spy,spz,ba,ca,pa,xc,yc,zc,ra,rb,phi,shift);
%% for Hamilton coordinates
xc_new = xc;
yc_new = yc;
ra_new = ra;
rb_new = rb;
phi_new = phi;
shift_new = shift;
Nitr = 5000;
ditr = 20;

%% drug parameters
gx = 0.4; % coordinate Hamilton parameter
gr = 20; % radius Hamilton parameter
%rr = 15; % radius restriction coefficient
gp = 0.13; % phi Hamilton parameter
gs = 0.6; % shift Hamilton parameter
constr = 0.1;
constrxy = 1;

n1 = num2str(kb);
n2 = num2str(kb_bott);
n3 = num2str(ks);
n4 = num2str(effective_r);
n5 = num2str(kz);
n6 = num2str(ka);
dir_name = 'Y:\tsygankov-lab\Nastya Zhurikhina\Simulations\Testing 3d new scaling\09.13.18\';
if type == 0
        s = mkdir(strcat(dir_name, 'wt_kb=', n1, ';kb_bott=', n2, ';ks=', n3, ';r=', n4, ';kz=', n5, ';ka=', n6)); %for Windows
        dir_name = strcat(dir_name, 'wt_kb=', n1, ';kb_bott=', n2, ';ks=', n3, ';r=', n4, ';kz=', n5, ';ka=', n6,'/');
elseif type == 1
        s = mkdir(strcat(dir_name, 'ccm1_kb=', n1, ';kb_bott=', n2, ';ks=', n3, ';r=', n4, ';kz=', n5, ';ka=', n6)); %for Windows
        dir_name = strcat(dir_name, 'ccm1_kb=', n1, ';kb_bott=', n2, ';ks=', n3, ';r=', n4, ';kz=', n5, ';ka=', n6,'/');
elseif type == 3
        s = mkdir(strcat(dir_name, 'ccm3_kb=', n1, ';kb_bott=', n2, ';ks=', n3, ';r=', n4, ';kz=', n5, ';ka=', n6)); %for Windows
        dir_name = strcat(dir_name, 'ccm3_kb=', n1, ';kb_bott=', n2, ';ks=', n3, ';r=', n4, ';kz=', n5, ';ka=', n6,'/');
end

test_info = [dir_name 'test_info.mat'];
save(test_info,'nc','np','nd','R','coord_len','V0','m','dl','dla','k','kl','kb','kb_bott','att','q','ks','ka','kz','t','Nitr','ditr','gx','gr','gp','gs','effective_r','effective_n');      

%% Iterations runing
for itr = 1:Nitr
     tic 
    %% generating prot lengths
    %[prot_len, es] = prot_lengths(np,nd,nc,k,kl,dl,prot_len,L,ba,ca,es,dla,ua,wa,ud,wd,mld,mla);
    %[prot_len,es] = prot_lengths_old(es,xc,yc,np,nd,nc,dl,prot_len,begx,begy,endx,endy,ca,ba,dla,R,effective_r,effective_n);
    [prot_len, es, ae] = prot_lengths(es,xc,yc,np,nd,nc,dl,prot_len,begx,begy,endx,endy,ca,ba,dla,R,effective_r,close_r,effective_n,close_n);
    %% checking for any attachments between cells or to the box
    [prot_len,endx,endy,endz,spx,spy,spz,ba,es,ca,pa] = prot_attachments(att,nc,np,nd,prot_len,begx,begy,begz,endx,endy,endz,spx,spy,spz,Box_coord,ba,es,ca,pa,ae,xc,yc,zc,ra,rb,phi,R,H,kb,kb_bott,dla);
    %check_length(itr) = prot_len(1);

    %% Hamilton gradient for all parameters
    parfor num = 1:nc
       %%{
        H1 = system_hamiltonian(1,num,nc,np,nd,prot_len,begx,begy,begz,endx,endy,endz,spx,spy,spz,ba,pa,ca,xc,yc,zc,ra,rb,phi,shift,V,N,ka,kz,ks,+q/2);
        H2 = system_hamiltonian(1,num,nc,np,nd,prot_len,begx,begy,begz,endx,endy,endz,spx,spy,spz,ba,pa,ca,xc,yc,zc,ra,rb,phi,shift,V,N,ka,kz,ks,-q/2);
        delta = t*(H1-H2)/(q*gx);
        xc_new(num) = xc(num)-delta*(abs(delta) < constrxy) - constrxy*sign(delta)*(abs(delta) >= constrxy); 

        H1 = system_hamiltonian(2,num,nc,np,nd,prot_len,begx,begy,begz,endx,endy,endz,spx,spy,spz,ba,pa,ca,xc,yc,zc,ra,rb,phi,shift,V,N,ka,kz,ks,+q/2);
        H2 = system_hamiltonian(2,num,nc,np,nd,prot_len,begx,begy,begz,endx,endy,endz,spx,spy,spz,ba,pa,ca,xc,yc,zc,ra,rb,phi,shift,V,N,ka,kz,ks,-q/2);
        delta = t*(H1-H2)/(q*gx);
        yc_new(num) = yc(num)-delta*(abs(delta) < constrxy) - constrxy*sign(delta)*(abs(delta) >= constrxy); 
        
        H1 = system_hamiltonian(3,num,nc,np,nd,prot_len,begx,begy,begz,endx,endy,endz,spx,spy,spz,ba,pa,ca,xc,yc,zc,ra,rb,phi,shift,V,N,ka,kz,ks,+q/2);
        H2 = system_hamiltonian(3,num,nc,np,nd,prot_len,begx,begy,begz,endx,endy,endz,spx,spy,spz,ba,pa,ca,xc,yc,zc,ra,rb,phi,shift,V,N,ka,kz,ks,-q/2);
        delta = t*(H1-H2)/(q*gr);
        ra_new(num) = ra(num)-delta*(abs(delta) < constr) - constr*sign(delta)*(abs(delta) >= constr); 
        
        H1 = system_hamiltonian(4,num,nc,np,nd,prot_len,begx,begy,begz,endx,endy,endz,spx,spy,spz,ba,pa,ca,xc,yc,zc,ra,rb,phi,shift,V,N,ka,kz,ks,+q/2);
        H2 = system_hamiltonian(4,num,nc,np,nd,prot_len,begx,begy,begz,endx,endy,endz,spx,spy,spz,ba,pa,ca,xc,yc,zc,ra,rb,phi,shift,V,N,ka,kz,ks,-q/2);
        delta = t*(H1-H2)/(q*gr);
        rb_new(num) = rb(num)-delta*(abs(delta) < constr) - constr*sign(delta)*(abs(delta) >= constr);   
        
        H1 = system_hamiltonian(5,num,nc,np,nd,prot_len,begx,begy,begz,endx,endy,endz,spx,spy,spz,ba,pa,ca,xc,yc,zc,ra,rb,phi,shift,V,N,ka,kz,ks,+q/2);
        H2 = system_hamiltonian(5,num,nc,np,nd,prot_len,begx,begy,begz,endx,endy,endz,spx,spy,spz,ba,pa,ca,xc,yc,zc,ra,rb,phi,shift,V,N,ka,kz,ks,-q/2);
        phi_new(num) = phi(num)-t*(H1-H2)/(q*gp);
        %disp([phi(num),phi_new(num)]);
        
        H1 = system_hamiltonian(6,num,nc,np,nd,prot_len,begx,begy,begz,endx,endy,endz,spx,spy,spz,ba,pa,ca,xc,yc,zc,ra,rb,phi,shift,V,N,ka,kz,ks,+q/2);
        H2 = system_hamiltonian(6,num,nc,np,nd,prot_len,begx,begy,begz,endx,endy,endz,spx,spy,spz,ba,pa,ca,xc,yc,zc,ra,rb,phi,shift,V,N,ka,kz,ks,-q/2);
        shift_new(num) = shift(num)-t*(H1-H2)/(q*gs); 
        %disp([shift(num),shift_new(num)]); 
        
    end
    
    xc = xc_new;
    yc = yc_new;
    ra = ra_new;
    rb = rb_new;
    phi = phi_new;
    shift = shift_new;
    [begx, begy, begz, endx, endy, endz, spx, spy, spz, L] = prot_coord(np,nc,nd,N,prot_len,Box_coord,spx,spy,spz,ba,ca,pa,xc,yc,zc,ra,rb,phi,shift);
    %% drawing cell  
%    %{
    if mod(itr,ditr) == 0 || itr == 1
        
        if mod(itr,ditr) == 0
            mat_num = itr/ditr + 1;
        end    
        if itr == 1
            mat_num = itr;
        end  
        %if itr == Nitr
        flnm=[dir_name num2str(mat_num) '.mat'];
        save(flnm,'R','H','np','nd','nc','begx','begy','begz','endx','endy','endz','spx','spy','spz','ca','ba','xc','yc','zc','ra','rb','phi','shift','V');       
        disp(V(1));
        %end

    end   
    toc
%   %}
end

return; 


end
