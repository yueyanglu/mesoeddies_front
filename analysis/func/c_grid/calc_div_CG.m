function divF = calc_div_CG(Fu,Fv,dyCu,dxCv,dxT,dyT,flg)
%
%   Calculate horizontal flux divergence on a C-grid. 
%   Fx shoule be on the U-grid, and Fy on the V-grid. 
%   The divergence is defined as the net outgoing transport across
%   the FOUR lateral boundaries divided by the area of the cell. 
%   divF ~ [FLUX/L].
%
%   Syntax: divF = calc_div_HYCOM(Fu,Fv,scuy,scvx,scp2)
%             
%         divF --    [nih-njh] p-grid          
%      Fu,dyCu --    [niu-nju] u-grid
%      Fv,dyCu --    [niv-njv] v-grid
% 

[niu,~] = size(Fu); 
[~,njv] = size(Fv); 
[nip,njp] = size(dxT);

dxyT2 = dxT .* dyT;

if niu ~= nip+1 || njv ~= njp+1
    error('Grid and flux dimensions do not match!');
end

[ic, jc] = deal(1:nip, 1:njp);

%%% enforce BC of no-norm flux. This only makes sense in closed domain! 
Fu([1 end],:) = 0;
Fv(:,[1 end]) = 0;

%%%
if flg == 1
    divF = (Fu(ic+1,jc).*dyCu(ic+1,jc) - Fu(ic,jc).*dyCu(ic,jc) + ...
        Fv(ic,jc+1).*dxCv(ic,jc+1) - Fv(ic,jc).*dxCv(ic,jc)) ./ dxyT2;
elseif flg == 2
    divF = (Fu(ic+1,jc)- Fu(ic,jc))./dxT + (Fv(ic,jc+1) - Fv(ic,jc))./dyT;
end
