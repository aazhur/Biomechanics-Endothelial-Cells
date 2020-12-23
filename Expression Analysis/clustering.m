X = readmatrix('clust_table.csv');
%ind = [1,3,4,5,7,8,9,11,12];
%{
M = M(ind,:);
M(1:3,:) = M(1:3,:)./(M(3,:)+10^(-6)); 
M(4:6,:) = M(4:6,:)./(M(6,:)+10^(-6)); 
M(7:9,:) = M(7:9,:)./(M(9,:)+10^(-6));
%}
R = [1:12]';
no_dims = 2;
init_dims = 12;
perplexity = 30;
mappedX = tsne(X, [], no_dims, init_dims, perplexity);
gscatter(mappedX(:,1), mappedX(:,2), R);

%R = 1:12;
%{
[ichi,schi] = fscchi2(M,R);
indi = find(schi(ichi));
ichi = ichi(indi);

[ittest,stest] = fsrftest(M,R);
indt = find(stest(ittest));
ittest = ittest(indt);

c = intersect(ichi,ittest);
%}
%{
[idx,scores] = fscmrmr(M,R);
ind = find(scores(idx));
idx = idx(ind);
%}
%{
idx = c;
scores = (stest(c) + schi(c))/2; 
[scores,i] = sort(scores,'descend');
idx = idx(i);

figure;
h = histogram(scores,5);
bins = h.BinEdges;
clust = discretize(scores,bins);
%}
%{
Groups = [idx', clust', scores'];
writematrix(Groups,'groups.csv');
%}
