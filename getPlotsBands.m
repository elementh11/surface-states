
%clearvars;

ef = 12.8320;
load bands_W90.mat;
b_W90 = dataw.ek - ef;
low_W90 = min(min(b_W90));
high_W90 = max(max(b_W90));
    
figure;

for jj = 1 : size(b_W90,2)
    band = b_W90(:,jj);
    hold all; plot(band,'-b','LineWidth',.5);
end

xlength = size(b_W90,1);
xlim([1 xlength]);
hold all; plot([1 xlength],[0 0],':k');

b_W90 = dataw.eks - ef;
for jj = 1 : size(b_W90,2)
    band = b_W90(:,jj);
    hold all; plot(band,'-r');
end
