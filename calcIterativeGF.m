function [G, GS, GSD] = calcIterativeGF2(H00, H01, omega)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%     USE TRANSFER MATRIX TO CALCULATE SELF-ENERGY
% J. Phys. F: Met. Phys 15(1985)851-858 M. P. Lopez Sancho, J. M. Lopez Sancho and J. Rubio
%     ------------------------------------------------------------
nterx = 100;
tol = 1.0e-16;
%     ------------------------------------------------------------
H00 = omega - H00;
H01 = omega - H01;
H10 = H01';
%     ------------------------------------------------------------
E = H00;
ES = E;
ESD = E;
    alpha = H01;
    beta = H10;
%     ------------------------------------------------------------
for m=1:nterx
    alphaS = alpha/E;
    betaS = beta/E;
    alphaSbeta = alphaS*beta;
    betaSalpha = betaS*alpha;
    E = E + alphaSbeta + betaSalpha;
    alpha = alphaS*alpha;
    beta = betaS*beta;
    ES = ES + alphaSbeta;
    ESD = ESD + betaSalpha;
    %     ------------------------------------------------------------
    conv = max(max(abs(alpha)));
    conv2 = max(max(abs(beta)));
    if(conv<tol)&&(conv2<tol); break,end
        %fprintf('t-matrix convergence %f %f %d\n', conv, conv2,m)

end % nterx
%     ------------------------------------------------------------
% if(conv>tol)||(conv2>tol)
%     error('bad t-matrix convergence %f %f', conv, conv2)
% end
%     ------------------------------------------------------------

G = inv(omega - E + 1i*tol);
GS = inv(omega - ES + 1i*tol);
GSD = inv(omega - ESD + 1i*tol);
