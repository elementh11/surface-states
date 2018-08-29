% ---- bulk DOS calculation
clearvars;

ns = 10;
ne = 1;
ef = 12.8320;
%elist = ef + linspace(-1, 1, ne);
elist = ef;
eta = 1e-15;
alat = 6.5706;
a1 = [   1.000000   0.000000   0.000000 ] * alat;
a2 = [   0.000000   1.797000   0.000000 ] * alat;
a3 = [   0.000000   0.000000   4.032000 ] * alat;

vol = abs(dot(cross(a1,a2),a3));
b1 = 2 * pi * cross(a2,a3) / vol;
b2 = 2 * pi * cross(a3,a1) / vol;
b3 = 2 * pi * cross(a1,a2) / vol;


load hfile.mat;
matrices = datah.matrices;
nrpts = datah.nrpts;
dim = datah.num_wann;

load kfile.mat;
kpoints = datak.kpoints;
nk = size(kpoints,1);
ldosb = zeros(nk, ne);
ldoss = zeros(nk, ne);
ldossd = zeros(nk, ne);



for kc = 1:nk
    k = kpoints(kc,1:3);
    realk = k(1)*b1 + k(2)*b2 + k(3)*b3;

    H00 = zeros(dim);  H01 = zeros(dim);

    for counter = 1:nrpts
        matrix = matrices(counter);
        delta = matrix.disp;
        if delta(3) == 0
            realdisp = delta(1) * a1 + delta(2) * a2 + delta(3) * a3;
            ham = matrix.ham;
            hcontr = (ham * exp(1i* sum(conj(realk).*realdisp))* matrix.deg);
            H00 = H00 + hcontr;
        end
        if delta(3) == 1
            realdisp = delta(1) * a1 + delta(2) * a2 + delta(3) * a3;
            ham = matrix.ham;
            hcontr = (ham * exp(1i* sum(conj(realk).*realdisp))* matrix.deg);
            H01 = H01 + hcontr;
        end
    end

    H00 = .5 * (H00 + H00');
    H01 = .5 * (H01 + H01');



    ec = 0;
    for e = elist
        ec = ec + 1;
        omega = e*eye(dim);
        [G, GS, GSD] = calcIterativeGF(H00, H01, omega);



        ldosb(kc, ec) = -imag(trace(G)) / pi;
        ldoss(kc, ec) = -imag(trace(GS)) / pi;
        ldossd(kc, ec) = -imag(trace(GSD)) / pi;




    end




end
dataw.ldosb = ldosb;
dataw.ldoss = ldoss;
dataw.ldossd = ldossd;
save('ldos_W90','dataw');
