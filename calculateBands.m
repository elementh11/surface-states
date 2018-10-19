% ---- bulk and surface bands calculation
clearvars;

ns=10;
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

load kpointslines.mat;
kpoints = datak.kpoints;
nk = size(kpoints,1);
energy = zeros(nk, dim);
energys = zeros(nk, dim * ns);



for kc = 1:nk
    k = kpoints(kc,1:3);
    realk = k(1)*b1 + k(2)*b2 + k(3)*b3;
    
    hamWk = zeros(dim);
    hams = zeros(dim * ns);
    for i=1:ns
        for j=1:ns
            hblock{i,j} = zeros(dim);
        end
    end
    
    for counter = 1:nrpts
        matrix = matrices(counter);
        delta = matrix.disp;
        realdisp = delta(1) * a1 + delta(2) * a2 + delta(3) * a3;
        ham = matrix.ham;
        hcontr = (ham * exp(1i* sum(conj(realk).*realdisp))* matrix.deg);
        hamWk = hamWk + hcontr;
        
        if abs(delta(3)) <= ns
            for i=1:ns
                for j=1:ns
                    if (i-j) == delta(3)
                        
                        hblock{i,j} = hblock{i,j} + hcontr;
                    end
                end
            end
        end
    end
    
    hamWk = .5 * (hamWk + hamWk');
    [~, ek] = eig(hamWk); ek = diag(ek);
    energy(kc,:) = ek;
    
    for i=1:ns
        for j=1:ns
            hams((dim*i-dim+1):(dim*i), (dim*j-dim+1):(dim*j)) = hblock{i,j};
        end
    end
    
    hams = .5 * (hams + hams');
    [~, eks] = eig(hams); eks = diag(eks);
    energys(kc,:) = eks;
    
    
end
dataw.ek = energy;
dataw.eks = energys;
save('bands_W90','dataw');


