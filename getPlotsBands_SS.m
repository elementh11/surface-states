
%clearvars;

ef = 12.8320;
load bands_W90.mat;
sz=3;
figure;

b_W90 = abs(dataw.eks) - ef;
xlength = size(b_W90,1);
for jj = 1 : size(b_W90,2)
    band = b_W90(:,jj);
    %hold all; plot(band,'-r');
    hold all; scatter(1:xlength,band,sz,'filled','r');
end

xlim([1 xlength]);
hold all; plot([1 xlength],[0 0],':k');

b_W90 = abs(dataw.ek) - ef;
for jj = 1 : size(b_W90,2)
    band = b_W90(:,jj);
    %hold all; plot(band,'-b','LineWidth',.5);
    hold all; scatter(1:xlength,band,sz,'filled','b');
end



b_W90 = abs(dataw.eksd) - ef;
for jj = 1 : size(b_W90,2)
    band = b_W90(:,jj);
    %hold all; plot(band,'-g');
    hold all; scatter(1:xlength,band,sz,'filled','g');
end
ylim([-1 1]);
