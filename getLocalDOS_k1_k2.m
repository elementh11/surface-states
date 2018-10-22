% ---- bulk and surface DOS calculation
clearvars;

ns = 7;%<---- max number of layers in the slab
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
ldosstop = zeros(nk, ne);
ldossbottom = zeros(nk, ne);



for kc = 1:nk
    k = kpoints(kc,1:3);
    realk = k(1)*b1 + k(2)*b2 + k(3)*b3;
    
    hamWk = zeros(dim);
    hamsurftop = zeros(dim * ns);
    hamsurfbottom = zeros(dim * ns);
    for i=1:ns
        for j=1:ns
            hblocktop{i,j} = zeros(dim);
            hblockbottom{i,j} = zeros(dim);
        end
    end
    
    for counter = 1:nrpts
        matrix = matrices(counter);
        delta = matrix.disp;
        realdisp = delta(1) * a1 + delta(2) * a2 + delta(3) * a3;
        ham = matrix.ham;
        hcontr = (ham * exp(1i* sum(conj(realk).*realdisp))* matrix.deg);
        hamWk = hamWk + hcontr;
        
        if delta(3) >= 0
            for i=(ns+1)/2:ns
                for j=(ns+1)/2:ns
                    if (i-j) == delta(3)
                        hblocktop{i,j} = hblocktop{i,j} + hcontr;
                    end
                end
            end
        else
            for i=1:(ns+1)/2
                for j=1:(ns+1)/2
                    if (i-j) == delta(3)
                        hblockbottom{i,j} = hblockbottom{i,j} + hcontr;
                    end
                end
            end
        end
    end
    
    hamWk = .5 * (hamWk + hamWk');
    for i=1:ns
        for j=1:ns
            hamsurftop((dim*i-dim+1):(dim*i), (dim*j-dim+1):(dim*j)) = hblocktop{i,j};
            hamsurfbottom((dim*i-dim+1):(dim*i), (dim*j-dim+1):(dim*j)) = hblockbottom{i,j};
        end
    end
    hamsurftop = .5 * (hamsurftop + hamsurftop');
    hamsurfbottom = .5 * (hamsurfbottom + hamsurfbottom');
    
    
    
    ec = 0;
    for e = elist
        ec = ec + 1;
        gb = inv(e*eye(dim) - hamWk + 1i*eta);
        ldosb(kc, ec) = -imag(trace(gb)) / pi;
        gstop = inv(e*eye(dim*ns) - hamsurftop + 1i*eta);
        ldosstop(kc, ec) = -imag(trace(gstop)) / pi;
        gsbottom = inv(e*eye(dim*ns) - hamsurfbottom + 1i*eta);
        ldossbottom(kc, ec) = -imag(trace(gsbottom)) / pi;
    end
    
    
    
    
end
dataw.ldosbulk = ldosb;
dataw.ldosstop = ldosstop;
dataw.ldossbottom = ldossbottom;
save('ldos_k1_k2','dataw');


