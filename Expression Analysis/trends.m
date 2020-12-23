T = csvread('data_values.csv');

C1 = [T(:,2)-T(:,1), T(:,3)-T(:,2), T(:,5)-T(:,4), T(:,6)-T(:,5)];
C2 = [T(:,1)-T(:,7), T(:,2)-T(:,8), T(:,3)-T(:,9), T(:,4)-T(:,7), T(:,5)-T(:,8), T(:,6)-T(:,9)];

C3 = sign(C2);
C3 = [sum(C3(:,1:3),2), sum(C3(:,4:6),2)];
C3 = (C3(:,1) == 3).*(C3(:,2) == -3) + (C3(:,1) == -3).*(C3(:,2) == 3);

C4 = sign(C1);
C4 = [sum(C4(:,1:2),2), sum(C4(:,3:4),2)];
C4 = (C4(:,1) < 0).*(C4(:,2) > 0) + (C4(:,1) > 0).*(C4(:,2) < 0);

Min = min(T(:,7:9),[],2); 
Max = max(T(:,7:9),[],2); 
C51 = sign([T(:,1)-Min, T(:,2)-Min, T(:,3)-Min, T(:,4)-Max, T(:,5)-Max, T(:,6)-Max]);
C52 = sign([T(:,1)-Max, T(:,2)-Max, T(:,3)-Max, T(:,4)-Min, T(:,5)-Min, T(:,6)-Min]);
C51 = [sum(C51(:,1:3),2), sum(C51(:,4:6),2)];
C52 = [sum(C52(:,1:3),2), sum(C52(:,4:6),2)];
C5 = C51+C52;
C5 = (C5(:,1) == 6).*(C5(:,2) == -6) + (C5(:,1) == -6).*(C5(:,2) == 6);

CB = find(C3.*C4);
 %CB = find(C5);

T = T(CB,:);
C6 = [T(:,1)./T(:,7), T(:,2)./T(:,8), T(:,3)./T(:,9), T(:,4)./T(:,7), T(:,5)./T(:,8), T(:,6)./T(:,9)];
C6 = (abs(100*(C6-1)) >= 25);
C6 = [sum(C6(:,1:3),2), sum(C6(:,4:6),2)];
C6 = (C6(:,1) == 3)+(C6(:,2) == 3);
C = find(C6);

IND = CB(C);
f = fopen('GeneIND1.txt','wt');
for i = 1:size(IND,1)
    fprintf(f,'%g\t',IND(i,:));
    fprintf(f,'\n');
end
fclose(f);

M = T(C,:);

r = 1:size(M,1);

%genes = char('ACO2','NA','C6orf48','NA','COL1A2','CRKL','CUX1','DAG1','EMC7','ERO1L','NA','LINC00116','MPP5','NEK1','PRKCSH','NA','UBA3');

figure;
for x = 1:size(M,1)
    subplot(10,10,x);
    hold on;
    scatter(1, M(x,1), 'r'); scatter(1, M(x,2), 'r'); scatter(1, M(x,3), 'r');
    scatter(1, M(x,4), 'g'); scatter(1, M(x,5), 'g'); scatter(1, M(x,6), 'g');
    scatter(1, M(x,7), 'b'); scatter(1, M(x,8), 'b'); scatter(1, M(x,9), 'b');
    %title(genes(x,:));
end    
