function [phi, psi, udiv, vdiv, urot, vrot] = helmholtz_decomp(u, v, ifwet, dxCu, dxCv, dyCu, dyCv, WRAP)
% 
% Helmholtz decomposition based on the function 'helmholtz_Cgrid'.
% Boundary points are dealt with in this function.
% 
% Input: all on the symmetric C-grid
%   u    [niu,nju]:
%   v    [niv,njv]:
%  ifwater [nih,njh]: true where (i,j) is inside the domain.
%   WRAP [logical 2 element vector]: whether to wrap in each dimension
%    
% Output:       
%     phi  p-grid
%    psiH  p-grid
%     psi  q-grid
% 
%     [urot_psi, vrot_psi] = psi2rot(psi, grid.dyCu, grid.dxCv);
%     [udiv_phi, vdiv_phi] = calc_GxGy_CG(phi,ones(size(phi)),grid.dxCu,grid.dyCu,grid.dxCv,grid.dyCv,0);
% 

% ---- check size
[nih,njh] = size(ifwet);
[niu,nju] = size(u);
[niv,njv] = size(v);
if nih ~= niu-1 || njh ~= njv-1 
    error('Fields must be defined on C-grid with [x,y] !!!');
end

% ----- NaN to zero
u(isnan(u)) = 0;
v(isnan(v)) = 0;

% ----- delete the N&E points (consistent with "helmholtz_Cgrid"), all onto p-size
u_p = u(1:end-1,:);
v_p = v(:,1:end-1);
e1u = dxCu(1:end-1,:);
e1v = dxCv(:,1:end-1);
e2u = dyCu(1:end-1,:);
e2v = dyCv(:,1:end-1);
% --- WARNING: the following (deleting S&W points) causes large errors!
% u_p = u(2:end,:);
% v_p = v(:,2:end);
% e1u = dxCu(2:end,:);
% e1v = dxCv(:,2:end);
% e2u = dyCu(2:end,:);
% e2v = dyCv(:,2:end);

% ----- decompose
[phi, psiH, psi, urot, vrot, udiv, vdiv] = helmholtz_Cgrid(u_p, v_p, ifwet, ...
    e1u, e1v, e2u, e2v, WRAP, true);

% ----- add N&E boundary points; thus, u is [niu,nju], v...
urot = cat(1, urot, zeros(1,nju));
vrot = cat(2, vrot, zeros(niv,1));
udiv = cat(1, udiv, zeros(1,nju));
vdiv = cat(2, vdiv, zeros(niv,1));

% psi on q-
if WRAP(1) && WRAP(2)
    psi = cat( 1, psiH, zeros(1,njh) );
    psi = cat( 2, psi, zeros(nih+1,1) );
end

% "helmholtz_Cgrid" solves (u,v) = grad[chi] - z X grad[psi]
% transfer it to (u,v) = grad[chi] + z X grad[psi]
psi = -psi;

end