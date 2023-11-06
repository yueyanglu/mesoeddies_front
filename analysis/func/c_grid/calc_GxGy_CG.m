function [Gx,Gy] = calc_GxGy_CG(c,h,dxCu,dyCu,dxCv,dyCv,flg)
%
% Calculate grid-weighted tracer gradients.
%   Gx [c or c*m] can be multiplied by the open face width, consistent with
%   the MOM6. If flg=0, not multiply, Gx~[c]; if 1, Gx~[face/dx*h*c]~[c*m] (MOM6).
% 
% Outputs:
%   Gx,Gy -- gradients on u/v points. The unit is controled by 'flg' 
% 
% Inputs (all defined on c-grid):
%   c -- tracer fld
%   h -- layer thickness [m]
%   flg -- flag for whether multiply gradient by grid len. '0' is NO [c],
%          others are YES [m*c].
% 
% Syntax: [Gx,Gy] = calc_GxGy_CG(c,h,dxCu,dyCu,dxCv,dyCv,flg)
% 

% -------------------------------------------- params
onemu = 1.e-12; % thickness of vanishing layer
[nih,njh] = size(c);
[niu,nju] = size(dxCu);
[niv,njv] = size(dxCv);
if nih ~= niu-1 || njh ~= njv-1 
    error('Fields must be defined on C-grid with [x,y] !!!');
end

% if flg=0, grad is not multiplied by open face width: dyCu and dxCv
if flg == 0
    dyCu = ones(size(dyCu));
    dxCv = ones(size(dxCv));
end

% -------------------------------------------- zonal grad, W&E are NaNs [u-grid]
% layer thickness at the interface by harmonic mean [u-grid]
h_harm_u = harmon( max(h(1:end-1,:),onemu), max(h(2:end,:),onemu) ); 
h_harm_u = cat(1, zeros(1,nju), h_harm_u, zeros(1,nju));

% delta_c along x-axis [u-grid]
dc_u = cat(1, zeros(1,nju), c(2:end,:)-c(1:end-1,:), zeros(1,nju)); 

% x tracer gradient, face/dx * h * delta_c
Gx = dyCu./dxCu .* h_harm_u .* dc_u; 

% -------------------------------------------- meri grad, S&N are NaNs [v-grid]
% layer thickness at the interface by harmonic mean [v-grid]
h_harm_v = harmon( max(h(:,1:end-1),onemu), max(h(:,2:end),onemu) ); 
h_harm_v = cat(2, zeros(niv,1), h_harm_v, zeros(niv,1));

% delta_c along y-axis [v-grid]
dc_v = cat(2, zeros(niv,1), c(:,2:end)-c(:,1:end-1), zeros(niv,1));

% y tracer gradient, face/dx * h * delta_c
Gy = dxCv./dyCv .* h_harm_v .* dc_v; 

end

% harmonic mean of layer thickiness at the interface
function h = harmon(x,y)
h = 2 * x .* y ./ (x + y);
end
