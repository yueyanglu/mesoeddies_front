function [means,eddys] = tracflx_MOM(c,us,ue,vs,ve,scp2,len,flxmethd)
% 
% 'tracflux_MOM.m' calculates the 2D large-scale tracer and tracer fluxes 
% fields from the given tracer field and fluxes. Functions needed:
%  smooth_geom_CG.m & calc_TFluxes_CG.m
% 
%   Input:
%      c    [nih,njh]: instantaneous tracer field
%      u    [niu,nju]: Zonal flx 
%      v    [niv,njv]: Meridional flx
%      len : window length of spatial running mean
%      sc*2: area of different grid cells [m2]
%      methd: method for calculating flux
%
%   Output:
%      cs : large-scale tracer
% 


%% check size of input

% must be 2D
if numel(size(c)) > 2
    error('c2flx function accepts 2D fields only !!');
end

% c,uflx,vflx must be in the same size
if size(c,1) ~= size(us,1)-1 || size(c,2) ~= size(vs,2)-1 
    error('Fields must be defined on C-grid with [x,y] !!!');
end

%%

%----------------------------------------------- <c> and c'
cs = smooth_geom_CG(c, scp2, len, len);
ce = c - cs;
    
%----------------------------------------------- tracer fluxes

% <u><c>
[uscs, vscs] = calc_TFluxes_CG(us, vs, cs, flxmethd);

% u'<c>
[uecs, vecs] = calc_TFluxes_CG(ue, ve, cs, flxmethd);

% <u>c'
[usce, vsce] = calc_TFluxes_CG(us, vs, ce, flxmethd);
    
% u'c'
[uece, vece] = calc_TFluxes_CG(ue, ve, ce, flxmethd);

%----------------------------------------------- assign to sstruct
means.cs = cs; 
means.uscs = uscs; means.vscs = vscs; 
% 
eddys.ce = ce; 
eddys.uecs = uecs; eddys.vecs = vecs; 
eddys.usce = usce; eddys.vsce = vsce; 
eddys.uece = uece; eddys.vece = vece; 

