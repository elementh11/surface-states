% ---- OPTICAL CONDUCTIVITY WITH KUBO FORMULA
%
% ------  initialize A for each frequency
%         skip calculation for \Delta E < e_thresh
%         calculates both bulk and surface states contributions
%
clearvars;
ef = 12.8620;
e_thresh = 1e-4;
e_thresh2 = 1e-3;
T = 0;
kb = 1.38e-23;  qelec = 1.6e-19; kt = kb / qelec * T;
delta_ = 5*1e-3;%--- need to figure it out the width of the broadening

freqlist = linspace(0, 0.25, 100);
dimf = length(freqlist);
ns = 3;
alat = 6.5706;
a1 = [   1.000000   0.000000   0.000000 ] * alat;  
a2 = [   0.000000   1.797000   0.000000 ] * alat;
a3 = [   0.000000   0.000000   4.032000 ] * alat;  

vol = abs(dot(cross(a1,a2),a3));
b1 = 2 * pi * cross(a2,a3) / vol;
b2 = 2 * pi * cross(a3,a1) / vol;
b3 = 2 * pi * cross(a1,a2) / vol;
prefac = norm(b1) * norm(b2) *norm(b3) / (2*pi)^3;
prefacs = norm(b1) * norm(b2) / (2*pi)^2;


load hfile.mat;
matrices = datah.matrices;
nrpts = datah.nrpts;
dim = datah.num_wann;


load kfile0.mat;%<----- change kfile name accordingly !!!
kpoints = datak.kpoints;
%weight = datak.weightlist;
nk = size(kpoints,1);

for dir1 = {'x','y'}
    for dir2 = {'x','y'}
        if (dir1{1}=='y' && dir2{1}=='x')
            continue
        else
            sigma_inter_kw.(dir1{1}).(dir2{1}) = zeros(nk, dimf);
            sigmas_inter_kw.(dir1{1}).(dir2{1}) = zeros(nk, dimf);
            sigma_intra_kw.(dir1{1}).(dir2{1}) = zeros(nk, dimf);
            sigmas_intra_kw.(dir1{1}).(dir2{1}) = zeros(nk, dimf);
        end
    end
end
A = zeros(dim);

for kc = 1:nk
    k = kpoints(kc,1:3);
    realk = k(1)*b1 + k(2)*b2 + k(3)*b3;

    hamWk = zeros(dim);
    hams = zeros(dim * ns);
    for i=1:ns
        for j=1:ns
            hblock{i,j} = zeros(dim);
            for dir1 = {'x','y'}
                vblock.(dir1{1}){i,j} = zeros(dim);
            end
        end
    end
    
    for dir1 = {'x','y'}
        velocity.(dir1{1}) = zeros(dim);
        velocitys.(dir1{1}) = zeros(dim * ns);
    end
    
    for counter = 1:nrpts
        matrix = matrices(counter);
        delta = matrix.disp;
        realdisp = delta(1) * a1 + delta(2) * a2 + delta(3) * a3;
        ham = matrix.ham;
        hcontr = (ham * exp(1i* sum(conj(realk).*realdisp))* matrix.deg);
        vxcontr = 1i * realdisp(1) * hcontr;
        vycontr = 1i * realdisp(2) * hcontr;
        hamWk = hamWk + hcontr;
        velocity.x = velocity.x + vxcontr;
        velocity.y = velocity.y + vycontr;

        if abs(delta(3)) <= ns
            for i=1:ns
                for j=1:ns
                    if (i-j) == delta(3)
                        hblock{i,j} = hblock{i,j} + hcontr;
                        vblock.x{i,j} = vblock.x{i,j} + vxcontr;
                        vblock.y{i,j} = vblock.y{i,j} + vycontr;
                    end
                end
            end
        end
    end
    
    hamWk = .5 * (hamWk + hamWk');
    [Uk, ek] = eig(hamWk); ek = diag(ek);
    for dir1 = {'x','y'}
        velocity.(dir1{1}) = Uk' * velocity.(dir1{1}) * Uk;
        velocity.(dir1{1}) = .5 * (velocity.(dir1{1}) + velocity.(dir1{1})');
    end
    
    for i=1:ns
        for j=1:ns
            hams((dim*i-dim+1):(dim*i), (dim*j-dim+1):(dim*j)) = hblock{i,j};
            for dir1 = {'x','y'}
                velocitys.(dir1{1})((dim*i-dim+1):(dim*i), (dim*j-dim+1):(dim*j)) = vblock.(dir1{1}){i,j};
            end
        end
    end
    hams = .5 * (hams + hams');
    [Uks, eks] = eig(hams); eks = diag(eks);
    for dir1 = {'x','y'}
        velocitys.(dir1{1}) = Uks' * velocitys.(dir1{1}) * Uks;
        velocitys.(dir1{1}) = .5 * (velocitys.(dir1{1}) + velocitys.(dir1{1})');
    end
    
    
    
    B = ((abs(ek - ef) < e_thresh2) + 0);
    Bs = ((abs(eks - ef) < e_thresh2) + 0);
    for dir1 = {'x','y'}
        for dir2 = {'x','y'}
            if (dir1{1}=='y' && dir2{1}=='x')
                continue
            else
                wplasma_k.(dir1{1}).(dir2{1})(kc) =...
                    sum(diag(velocity.(dir1{1})).*B.*diag(velocity.(dir2{1})));...
                    %* weight(kc);
                wplasmas_k.(dir1{1}).(dir2{1})(kc) =...
                    sum(diag(velocitys.(dir1{1})).*Bs.*diag(velocitys.(dir2{1})));
            end
        end
    end
    
    fc = 0;
    for freq = freqlist
        fc = fc + 1;
        A = zeros(dim);
        for n = 1:dim
            if ek(n) > ef
                continue
            else
                %fnk = 1 / (exp((ek(n) - ef)/kt) + 1);
                for m = 1:dim
                    if (n == m) || (ek(m) <= ef) || (abs(ek(n)-ek(m)) < e_thresh)
                        continue
                    else
                        %fmk = 1 / (exp((ek(m) - ef)/kt) + 1);
                        A(n,m) = 1 / ((ek(n)-ek(m))*(ek(n)-ek(m) + freq + 1i*delta_));
                            
                    end
                end
            end
        end
        
        As = zeros(dim * ns);
        for n = 1:dim*ns
            if eks(n) > ef
                continue
            else
                %fnk = 1 / (exp((ek(n) - ef)/kt) + 1);
                for m = 1:dim*ns
                    if (n == m) || (eks(m) <= ef) || (abs(eks(n)-eks(m)) < e_thresh)
                        continue
                    else
                        %fmk = 1 / (exp((ek(m) - ef)/kt) + 1);
                        As(n,m) = 1 / ((eks(n)-eks(m))*(eks(n)-eks(m) + freq + 1i*delta_));
                            
                    end
                end
            end
        end
        
        for dir1 = {'x','y'}
            for dir2 = {'x','y'}
                if (dir1{1}=='y' && dir2{1}=='x')
                    continue
                else
                    sigma_inter_kw.(dir1{1}).(dir2{1})(kc, fc) =...
                        sum(sum(velocity.(dir1{1}).*A.*velocity.(dir2{1})'));...
                        %* weight(kc);
                    sigmas_inter_kw.(dir1{1}).(dir2{1})(kc, fc) =...
                        sum(sum(velocitys.(dir1{1}).*As.*velocitys.(dir2{1})'));...
                    sigma_intra_kw.(dir1{1}).(dir2{1})(kc, fc) =...
                        wplasma_k.(dir1{1}).(dir2{1})(kc) / (freq + 1i*delta_);
                    sigmas_intra_kw.(dir1{1}).(dir2{1})(kc, fc) =...
                        wplasmas_k.(dir1{1}).(dir2{1})(kc) / (freq + 1i*delta_);
                end
            end
        end
        
    end
    %     if mod(kc, 5000)==0
    %         prog = fopen('progress.txt','a');
    %         fprintf(prog,'calculating for k=%d/%d \n',kc,nk);
    %         fclose(prog);
    %     end
    fprintf('calculation for k = %d / %d ... done \n',kc,nk);
end

for dir1 = {'x','y'}
    for dir2 = {'x','y'}
        if (dir1{1}=='y' && dir2{1}=='x')
            continue
        else
            data.sigma_inter.(dir1{1}).(dir2{1}) = (sigma_inter_kw.(dir1{1}).(dir2{1}))*1i*prefac;
            data.sigma_intra.(dir1{1}).(dir2{1}) = (sigma_intra_kw.(dir1{1}).(dir2{1}))*1i*prefac;
            data.wplasma.(dir1{1}).(dir2{1}) = (wplasma_k.(dir1{1}).(dir2{1}));
            data.sigmas_inter.(dir1{1}).(dir2{1}) = (sigmas_inter_kw.(dir1{1}).(dir2{1}))*1i*prefacs;
            data.sigmas_intra.(dir1{1}).(dir2{1}) = (sigmas_intra_kw.(dir1{1}).(dir2{1}))*1i*prefacs;
            data.wplasmas.(dir1{1}).(dir2{1}) = (wplasmas_k.(dir1{1}).(dir2{1}));
        end
    end
end

save('opt_cond_SS','data','-v7.3');%<---- save with appropriate filename !!!


