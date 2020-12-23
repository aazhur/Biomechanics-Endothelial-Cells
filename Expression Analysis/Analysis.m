file = tdfread('genes.txt');
%{
ccm1 = file.ccm1_fpkm./file.NT_fpkm;
ccm3 = file.ccm3_fpkm./file.NT_fpkm;

ids1 = find((ccm1 >= 2).*(ccm3 >= 2) + (ccm1 <= 0.5).*(ccm3 <= 0.5));
ids2 = find((ccm1 >= 2).*(ccm3 <= 0.5) + (ccm1 <= 0.5).*(ccm3 >= 2));

names1 = file.gene_name(ids1,:);
names2 = file.gene_name(ids2,:);

dlmwrite('EID2.txt',names2,'delimiter','');
dlmwrite('EID1.txt',names1,'delimiter','');
%}

wnt = tdfread('Smad.txt');
ind = find(ismember(file.gene_name, wnt.Name, 'rows'));
wnt.ccm1 = file.ccm1_fpkm(ind);
wnt.ccm3 = file.ccm3_fpkm(ind);
wnt.wt = file.NT_fpkm(ind);

figure;
grid minor;
axis tight;
hold on;
N = length(wnt.Name);
scatter(linspace(1,N,N),wnt.ccm1./wnt.wt, 'r');
scatter(linspace(1,N,N),wnt.ccm3./wnt.wt, 'b');
scatter(linspace(1,N,N),wnt.wt./wnt.wt, 'g');

