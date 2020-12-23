function cell_simulation_loop(itr1,itr2)

%% variables and parameters used
% nc - number of cells
% xc - cell centers' x_coord
% yc - cell centers' y_coord
% ra - main elliptic diameter
% S0 - initial mean area
% S - area of cell
% phi - angel of cell rotation
% r0 - initial ellipse radius
% Xo - matrix of cells' exterior x-coordinates
% Yo - matrix of cells' exterior x-coordinates

% np - number of protrusions
% prot_len - lengths and attachment status of every protrusion
% ca - cell attachments of protrusions
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

%% box
    
    edge_dist = 50;
    coord_len = 650 + edge_dist;  

    xlim = [-coord_len, +coord_len];
    ylim = [-coord_len, +coord_len];
    s = subplot(1,1,1);
    set(s,'Box','on');
    set(s,'XLim',xlim,'YLim',ylim);
    set(gca,'Color', [0.3,0.3,0.3]);

    teta = 0:360;
    R = coord_len - edge_dist;
    disp(R);
    x_box = R*cosd(teta); y_box = R*sind(teta);

    %% main parameters of ellipse
    np = 16;
    nc = 1600;
    S0 = 40*pi;
    S = S0*ones(1 , nc);
    ra = sqrt(S0/pi)*ones(1 , nc);
    phi = 90*rand(1 , nc);
    shift = zeros(1, nc);
    r0 = sqrt(S(1)/pi);
    dl = 0.42;
    m = 1.2*r0; %minimal initial cell distance
    [xc,yc] = random_initial_placing(nc,r0,m,R);

    %% Hamilton gradient
    q = 0.01;
    ka = 5.0;
    ks = 0.5;
    t = 0.025;

    %% protrusions length parameters
    ua = 0.1;   % (mean-length) rate of switching to "growth" mode (when attached)
    wa = 0.9;   % (mean-length) rate of switching to "shrink" mode (when attached)
    
    ud = 0.5;   % (mean-length) rate of switching to "growth" mode (when dettached)
    wd = 0.5;   % (mean-length) rate of switching to "shrink" mode (when dettached)
    
    mld = 12;  %mean length when detached; 
    mla = 0;  %mean length when attached

    %% portrusions main parameters
    k = 0.1;
    kl = 3.5;
    l0 = 0;
    dla = 0.21; % dl when attached
    % set initial protrusions length and attachments
    prot_len = l0*ones(nc , np);

    %% springs
    kb_ = exp(linspace(4.8,8.1,10));
    kb = kb_(10);
    att = 1; 
    spx = zeros(nc , np);
    spy = zeros(nc , np);
    ca = zeros(nc , np);
    ba = zeros(nc , np);
    pa = zeros(nc , np);
    es = ones(nc , np);

    [begx, begy, endx, endy, spx, spy, L] = prot_coord(np,nc,prot_len,spx,spy,ba,ca,pa,xc,yc,ra,phi,shift,S);

    %% for Hamilton coordinates
    xc_new = xc;
    yc_new = yc;
    ra_new = ra;
    phi_new = phi;
    shift_new = shift;
    Nitr = 4000;
    ditr = 20;

    %% drug parameters
    %p = randperm(1600);
    %p = (p > 1600*0);
    gx_ = exp(linspace(-0.92,1.76,Nitr));
    %gx = gx_(itr2)*p + 0.4*(1-p);  
    %gx = gx_(itr2);
    gr = 20; % radius Hamilton parameter
    %rr = 15; % radius restriction coefficient
    gp = 0.13; % phi Hamilton parameter
    gs = 0.6; % shift Hamilton parameter

    number = int2str(itr2+(itr1-1)*10);
    s = mkdir(strcat('Y:\tsygankov-lab\Nastya Zhurikhina\Simulations\Testing CCM1 Descend\testing_11.28.2017_', number)); %for Windows
    disp(s);
    dir_name = strcat('Y:\tsygankov-lab\Nastya Zhurikhina\Simulations\Testing CCM1 Descend\testing_11.28.2017_', number,'\');
    test_info = [dir_name 'test_info.mat'];
    save(test_info,'nc','np','R','coord_len','S0','r0','m','dl','dla','mld','mla','k','kl','kb','att','q','ks','ka','t','Nitr','ditr','gx_','gr','gp','gs');     

%% Iterations runing
for itr = 1:Nitr
     tic 
    %% generating prot lengths
    [prot_len, es] = prot_lengths(np,nc,k,kl,dl,prot_len,L,ba,ca,es,dla,ua,wa,ud,wd,mld,mla);
    
    %% checking for any attachments between cells or to the box
    [prot_len,endx,endy,spx,spy,ba,es,ca,pa] = prot_attachments(att,nc,np,prot_len,begx,begy,endx,endy,spx,spy,x_box,y_box,ba,es,ca,pa,xc,yc,ra,phi,S,R,kb,dla);
   
    %% Hamilton gradient for all parameters
    gx = gx_(itr);
    parfor num = 1:nc
        H1 = system_hamiltonian(1,num,nc,np,begx,begy,endx,endy,spx,spy,ba,pa,ca,xc,yc,ra,phi,shift,S,ka,ks,+q/2);
        H2 = system_hamiltonian(1,num,nc,np,begx,begy,endx,endy,spx,spy,ba,pa,ca,xc,yc,ra,phi,shift,S,ka,ks,-q/2);
        xc_new(num) = xc(num)-t*(H1-H2)/(q*gx);
        %disp([xc(num),xc_new(num)]);
        %disp(t*(H1-H2)/(q*gx));
        
        H1 = system_hamiltonian(2,num,nc,np,begx,begy,endx,endy,spx,spy,ba,pa,ca,xc,yc,ra,phi,shift,S,ka,ks,+q/2);
        H2 = system_hamiltonian(2,num,nc,np,begx,begy,endx,endy,spx,spy,ba,pa,ca,xc,yc,ra,phi,shift,S,ka,ks,-q/2);
        yc_new(num) = yc(num)-t*(H1-H2)/(q*gx);
        %disp(t*(H1-H2)/(q*gx));
        
        H1 = system_hamiltonian(3,num,nc,np,begx,begy,endx,endy,spx,spy,ba,pa,ca,xc,yc,ra,phi,shift,S,ka,ks,+q/2);
        H2 = system_hamiltonian(3,num,nc,np,begx,begy,endx,endy,spx,spy,ba,pa,ca,xc,yc,ra,phi,shift,S,ka,ks,-q/2);
        ra_new(num) = ra(num)-t*(H1-H2)/(q*gr);
        %{
        if ra_new(num) > rr
            ra_new(num) = rr;
        end    
        if S0/(ra_new(num)*pi) > rr
            ra_new(num) = S0/(rr*pi);
        end  
       %}
        H1 = system_hamiltonian(4,num,nc,np,begx,begy,endx,endy,spx,spy,ba,pa,ca,xc,yc,ra,phi,shift,S,ka,ks,+q/2);
        H2 = system_hamiltonian(4,num,nc,np,begx,begy,endx,endy,spx,spy,ba,pa,ca,xc,yc,ra,phi,shift,S,ka,ks,-q/2);
        phi_new(num) = phi(num)-t*(H1-H2)/(q*gp);
        %disp([phi(num),phi_new(num)]);
        
        H1 = system_hamiltonian(5,num,nc,np,begx,begy,endx,endy,spx,spy,ba,pa,ca,xc,yc,ra,phi,shift,S,ka,ks,+q/2);
        H2 = system_hamiltonian(5,num,nc,np,begx,begy,endx,endy,spx,spy,ba,pa,ca,xc,yc,ra,phi,shift,S,ka,ks,-q/2);
        shift_new(num) = shift(num)-t*(H1-H2)/(q*gs); 
        %disp([shift(num),shift_new(num)]); 
        
    end
    xc = xc_new;
    yc = yc_new;
    ra = ra_new;
    phi = phi_new;
    shift = shift_new;
    
    [begx, begy, endx, endy, spx, spy,L] = prot_coord(np,nc,prot_len,spx,spy,ba,ca,pa,xc,yc,ra,phi,shift,S);
%{   
    cla;
    hold on;
    axis equal;

    draw_box(R, ba, np, nc);
    draw_prot_alt(np,nc,begx,begy,endx,endy,spx,spy,ca,ba,xc,yc,ra,phi,S);
    drawnow;
    pause(0.05);
 %}   
    
    %% drawing cell  
%    %{
    if mod(itr,ditr) == 0 || itr == 1
        
        if mod(itr,ditr) == 0
            mat_num = itr/ditr + 1;
        end    
        if itr == 1
            mat_num = itr;
        end  
        
        if itr == Nitr
            flnm=[dir_name num2str(mat_num) '.mat'];
            save(flnm,'R','es','np','nc','begx','begy','endx','endy','spx','spy','ca','ba','xc','yc','ra','phi','S','shift');        
        end
    end
        
    toc
%   %}
end
%draw_frames_Loop;
    
return; 
end
