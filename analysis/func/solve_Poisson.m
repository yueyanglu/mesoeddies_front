function [phi,dphidx,dphidy] = solve_Poisson(src, area, h, dxCu,dyCu,dxCv,dyCv,flg,ifdecompL)
% 
% Inverse del2{PHI} = div{del{PHI}} = src to obtain PHI; and returns del{PHI}
% PHI = L^(-1) * SRC * area
% 
% Recomanded: h = ones, flg = 0 
% 
if nargin < 9
    ifdecompL = 0;
end

[ni, nj] = size(src);

% Laplacian coef [dyCu * h / dxCu]
L = coef_Laplacian_CG(h,dxCu,dyCu,dxCv,dyCv); 
if ifdecompL
    dL = decomposition(L,'cod');
else
    dL = L;
end

% [src * m2]
F1d = reshape(src.*area, [ni*nj 1]);  

% phi [src * m2/L]
phi1d =  dL \ F1d; % pinv(dL) * F1d,  dL \ F1d or lsqminnorm(dL,F1d) - SLOW!
phi = reshape(phi1d, [ni, nj]);

% del{PHI}
[dphidx, dphidy] = calc_GxGy_CG(phi,h,dxCu,dyCu,dxCv,dyCv,flg);
% enforce no-normal-flux BC
dphidx([1 end],:) = 0;  dphidy(:,[1 end]) = 0;




