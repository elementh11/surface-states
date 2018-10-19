
%---- generates k-points along symmetry lines for Wannier bandstructure
%---- calculation

nk = 50;

% GY = [linspace(0.0, 0.0, nk);
%       linspace(0.0, 0.5, nk);
%       linspace(0.0, 0.0, nk)]';
% GY(nk, :) = [];
% 
% YS = [linspace(0.0, 0.5, nk);
%       linspace(0.5, 0.5, nk);
%       linspace(0.0, 0.0, nk)]';
% YS(nk, :) = [];
% 
% SX = [linspace(0.5, 0.5, nk);
%       linspace(0.5, 0.0, nk);
%       linspace(0.0, 0.0, nk)]';
% SX(nk, :) = [];

XG = [linspace(0.5, 0.0, nk);
      linspace(0.0, 0.0, nk);
      linspace(0.0, 0.0, nk)]';
%XG(nk, :) = [];

% GZ = [linspace(0.0, 0.0, nk);
%       linspace(0.0, 0.0, nk);
%       linspace(0.0, 0.5, nk)]';
% GZ(nk, :) = [];
% 
% ZU = [linspace(0.0, 0.0, nk);
%       linspace(0.0, 0.5, nk);
%       linspace(0.5, 0.5, nk)]';
% ZU(nk, :) = [];
% 
% UR = [linspace(0.0, 0.5, nk);
%       linspace(0.5, 0.5, nk);
%       linspace(0.5, 0.5, nk)]';
% UR(nk, :) = [];
% 
% RT = [linspace(0.5, 0.5, nk);
%       linspace(0.5, 0.0, nk);
%       linspace(0.5, 0.5, nk)]';
% RT(nk, :) = [];
% 
% TZ = [linspace(0.5, 0.0, nk);
%       linspace(0.0, 0.0, nk);
%       linspace(0.5, 0.5, nk)]';
% %TZ(nk, :) = [];



%klist = [GY; YS; SX; XG; GZ; ZU; UR; RT; TZ];
klist = XG;
datak.kpoints = klist;
save('kpointslines','datak');
















