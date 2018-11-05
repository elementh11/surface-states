% ---- (top / bottom) surface LDOS calculation using iterative Green function method
clearvars;

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
ab = zeros(nk, ne);
as = zeros(nk, ne);
asd = zeros(nk, ne);



for kc = 1:nk
    k = kpoints(kc,1:3);
    realk = k(1)*b1 + k(2)*b2 + k(3)*b3;
    %%%% build Hamiltonians of principal layers H00 and H01
    
    
    
    H00 = zeros(dim);
    H01 = zeros(dim);
    
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
    
    %%%%%%%% call iterative method to get G
    
    ec = 0;
    for energy = elist
        ec = ec + 1;
        omega = energy*eye(dim);
        [G, GS, GSD] = getIterativeGreenFunction(H00, H01, omega);

        %%%%%%% calculate spectral functions A(k, omega)

        ab(kc, ec) = -imag( trace (G)) / pi;
        as(kc, ec) = -imag( trace (GS)) / pi;
        asd(kc, ec) = -imag( trace (GSD)) / pi;
    end
    
    %%%%% insert progress indicator
    if mod(kc, 5000)==0
        prog = fopen('progress.txt','a');
        fprintf(prog,'calculating for k=%d/%d \n',kc,nk);
        fclose(prog);
    end
    
end

dataw.ldosb = ab;
dataw.ldoss = as;
dataw.ldossd = asd;
save('ldos_k1_k2','dataw');


