function [fu3d_sm,fv3d_sm] = smooth_flds(grid,win_len,fu3d,fv3d,ifdecomp)
%[fu3d_cs,fv3d_cs] = smooth_flds(grid,gridcs,fu3d,fv3d,ifdecomp)
% This function computes ONE snapshot (for 3d flds)!
% This is similar to "csgrain_flds.m" except that not projected onto a
% coarser grid.
% 
% Input: 
%    fu3d - u*h*L, will be scaled by dyCu before decomposition
%    ifdecomp - "0" for simple filtering; "1" for decomp-based
% 
% Spatially smoothing the fields (h, & flux F) and save them.
% This new SM method preserves the divergence of F, which is done by
% 
% 1. decompose F into DIV (phi) and ROT (psi) components: 
%      F = del{phi} + z X del{psi} = [dphi_dx, dphi_dy] + [-dpsi_dy, dpsi_dx]
% 
% 2. cs-graining div{F}, from which the DIV comp is obtianed by a Poisson eqn:
%      del2{phi_CS} = < div{F} >
% 
% 3. DIV comp of the cs-grained flux F_CS is defined as:
%      F_CS_DIV = del{phi_CS}
% 
% 4. ROT comp of F_CS is defined as: 
%      F_CS_ROT = z X del{ <psi> } = [-dpsi_dy, dpsi_dx]
% 
% 5. F_CS = F_CS_DIV + F_CS_ROT
% 

%% params
nk = size(fu3d,3);

[dxCv, dyCu] = deal(grid.dxCv, grid.dyCu);
[dxT, dyT] = deal(grid.dxT, grid.dyT); 
dxyT = grid.dxT .* grid.dyT;
%
ifpad = 1;
ifverb = 1;

%% do the snapshot

%----------------------- 3D vars to be obtained at each step
fu3d_sm = NaN * zeros(grid.niu, grid.nju, nk);
fv3d_sm = NaN * zeros(grid.niv, grid.njv, nk);

%----------------------- loop over k
for ik = 1:nk
    fprintf(1,'\nDo k = %02d of %02d\n', ik,nk);
    
    %---------- read flds
    % use uh = m2/s
    fu = fu3d(:,:,ik) ./ grid.dyCu;
    fv = fv3d(:,:,ik) ./ grid.dxCv;
    
    % div [fu/m] and curl [fu/m] of flux
    divF = calc_div_CG(fu,fv,dyCu,dxCv,dxT,dyT,1);

    %---------------------------------------------------
    % spatial smoothing
    %---------------------------------------------------
    [fu_sm,fv_sm] = smooth_uv(grid,win_len,fu,fv,ifdecomp,ifverb);

    %---- check the divergence [fu/m ~ m/s]
    divF_sm = smooth_geom_CG(divF,dxyT,win_len,win_len,2,ifpad); % div [1/s]
    div_F_sm = calc_div_CG(fu_sm,fv_sm,grid.dyCu,grid.dxCv,grid.dxT,grid.dyT,1);
    %     curl_F_cs = calc_curl_CG(fu_cs,fv_cs,grid_cs.dxBu,grid_cs.dyBu);
    fprintf( 'RMS(div_F_sm-src_cs) / RMS(src_cs) = %.1e / %.1e \n', ...
        rms(div_F_sm(:)-divF_sm(:)), rms(divF_sm(:),'omitnan') );
    
    %-----------------------------------------------------------
    %    assign 
    %-----------------------------------------------------------
    fu3d_sm(:,:,ik) = fu_sm .* grid.dyCu;
    fv3d_sm(:,:,ik) = fv_sm .* grid.dxCv;
    
end % k


%% plot
%{
fac = 1e6;

figure
% 
subplot(221)
hd = pcolor(fuD_phi_cs'); set(hd, 'EdgeColor', 'none');
cmocean('balance')
caxis([-.1 .1]*fac)
colorbar
title('CS grained uh')
% 
subplot(222)
hd = pcolor(fuR_psi_cs'); set(hd, 'EdgeColor', 'none');
cmocean('balance')
caxis([-1 1]*fac)
colorbar
title('<U>*: U to be used')
% 
subplot(223)
hd = pcolor(fu_cs'); set(hd, 'EdgeColor', 'none');
cmocean('balance')
caxis([-1 1]*fac)
colorbar
title('div comp of <U>*')
% 
subplot(224)
hd = pcolor(fu'); set(hd, 'EdgeColor', 'none');
cmocean('balance')
caxis([-1 1]*fac)
colorbar
title('rot comp of <U>*')

%%
figure
subplot(221)
hd = pcolor(src'); set(hd, 'EdgeColor', 'none');
cmocean('balance')
caxis([-5e-7 5e-7])
colorbar
title('<dh/dt>')
% 
subplot(222)
hd = pcolor(src_cs'); set(hd, 'EdgeColor', 'none');
cmocean('balance')
% caxis([-5e-4 5e-4])
caxis([-5e-7 5e-7])
colorbar
% 
subplot(223)
hd = pcolor(div_F_cs'); set(hd, 'EdgeColor', 'none');
cmocean('balance')
% caxis([-5e-4 5e-4])
caxis([-5e-7 5e-7])
colorbar
%}


end
