

%clearvars;

nk1 = 300;
nk2 = 300;
nk3 = 1;
load ldos_W90.mat;
doss = dataw.ldoss(:, 1);
dosb = dataw.ldosb(:, 1);
dossd = dataw.ldossd(:, 1);


kc = 0;
for k1counter = 1:nk1
    for k2counter = 1:nk2
        for k3counter = 1:nk3
            kc = kc+1;
            surf(k1counter,k2counter) = doss(kc);
            surfd(k1counter,k2counter) = dossd(kc);
            bulk(k1counter,k2counter) = dosb(kc);
        end
    end
end






