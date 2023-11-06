function [ut,vt] = calc_TFluxes_CG(u,v,t,method)
%   Calculates lateral advective tracer fluxes on C-grid
% 
%   SYNTAX: [ut,vt]=calc_TFluxes_HYCOM(u,v,t,method)
%           ut, vt -----  [niu/v-nju/v] [u*t]
%           u,v, --- [niu/v-njhu/v] velocities
%           t --- [nih-njh] tracer
%           method --- '2nd_order' (default) or '4th_order'

if nargin==3
    method='2nd_order';
    disp('Default 2nd order advection is used!')
end

% initial
[nih,njh] = size(t);
[niu,~] = size(u); 
[~,njv] = size(v); 
if niu ~= nih+1 || njv ~= njh+1
    error('Flux and t dimensions do not match!');
end

ut = NaN * zeros(size(u));
vt = NaN * zeros(size(v));

% --------- 
if strcmpi(method,'2nd_order')
%
%   Calculate 2nd order tracer fluxes
%   
    ut(2:end-1,:) = .5*u(2:end-1,:) .* (t(2:end,:) + t(1:end-1,:));
    vt(:,2:end-1) = .5*v(:,2:end-1) .* (t(:,2:end) + t(:,1:end-1));

elseif strcmpi(method,'4th_order')
%
%   Calculate 4th order tracer fluxes
%    Needs 4 h's and 1 u: 2 h's on the left and 2 on the right of the u!
%    x--o--x--o--x--o--x--o--x
%       h     h  u  h     h
% 
    iB = 3:niu-2;
    ic = 3:nih-1;
    ut(iB,:) = .5*u(iB,:) .* ( 1.125*t(ic,:)+1.125*t(ic-1,:)-0.125*t(ic+1,:)-0.125*t(ic-2,:) );

    jB = 3:njv-2;
    jc = 3:njh-1;
    vt(:,jB) = .5*v(:,jB) .* ( 1.125*t(:,jc)+1.125*t(:,jc-1)-0.125*t(:,jc+1)-0.125*t(:,jc-2) );

    %--- 2nd boundary points use 2nd order method
    ut([2 end-1],:) = .5*u([2 end-1],:) .* (t([1 end-1],:) + t([2 end],:));
    vt(:,[2 end-1]) = .5*v(:,[2 end-1]) .* (t(:,[1 end-1]) + t(:,[2 end]));

end

%--- 1st boundary points use an up-wind scheme, simply u*h
% can be ignored since u is NaN at the closed boundary
ut([1 end],:) = u([1 end],:) .* t([1 end],:);
vt(:,[1 end]) = v(:,[1 end]) .* t(:,[1 end]);

