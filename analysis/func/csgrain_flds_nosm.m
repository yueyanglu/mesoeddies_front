function [fu3d_cs,fv3d_cs] = csgrain_flds_nosm(grid,grid_cs,fu3d,fv3d,ifdecomp)
%[fu3d_cs,fv3d_cs,h3d_cs] = csgrain_flds_func(grid,gridcs,fu3d,fv3d,h3d)
% This function computes ONE snapshot (for 3d flds)!
% 
% Input: 
%    fu3d - u*h*L, will be scaled by dyCu before decomposition
%    ifdecomp - "0" for simple filtering u/v component; "1" for decomp-based

% % Coarse-grain the fields (h, & flux F) and save them.
% This new CS method preserves the divergence of F, which is done by
% 
% 1. decompose F into DIV (phi) and ROT (psi) components: 
%      F = del{phi} + z X del{psi} = [dphi_dx, dphi_dy] + [-dpsi_dy, dpsi_dx]
% 
% 2. cs-graining (sub-sample ONLY) div{F}, from which the DIV comp is obtianed by a Poisson eqn:
%      del2{phi_CS} = < div{F} >
% 
% 3. DIV comp of the cs-grained (sub-sampled) flux F_CS is defined as:
%      F_CS_DIV = del{phi_CS}
% 
% 4. ROT comp of F_CS is defined as: 
%      F_CS_ROT = z X del{ <psi> } = [-dpsi_dy, dpsi_dx]
% 
% 5. F_CS = F_CS_DIV + F_CS_ROT
% 

%% params
nk = size(fu3d,3);
win_len = grid.nih / grid_cs.nih;
if mod(win_len,1)~=0;  error('grid and grid_cs not scaled well!!!'); end


onesH = ones(grid.nih,grid.njh);
[dxCu, dxCv, dyCu, dyCv] = deal(grid.dxCu, grid.dxCv, grid.dyCu, grid.dyCv);
[dxBu, dyBu, dxT, dyT] = deal(grid.dxBu, grid.dyBu, grid.dxT, grid.dyT);    
% 
onesHcs = ones(grid_cs.nih,grid_cs.njh);
areaTcs = grid_cs.dxT .* grid_cs.dyT;

%% do the snapshot

%----------------------- 3D vars to be obtained at each step
fu3d_cs = NaN * zeros(grid_cs.niu, grid_cs.nju, nk);
fv3d_cs = NaN * zeros(grid_cs.niv, grid_cs.njv, nk);

%----------------------- loop over k
for ik = 1:nk
    fprintf(1,'\nDo k = %02d of %02d\n', ik,nk);
    
    %---------- read flds
    % use uh = m2/s
    fu = fu3d(:,:,ik) ./ dyCu;
    fv = fv3d(:,:,ik) ./ dxCv;
    
    % div [fu/m] and curl [fu/m] of flux
    divF = calc_div_CG(fu,fv,dyCu,dxCv,dxT,dyT,1);
    curl = calc_curl_CG(fu,fv,dxBu,dyBu);
    
    if ifdecomp % decompose-based filtering
        %
        %         curl_cs = sepblockfun(curl(1:end-1,1:end-1),[win_len,win_len],'nanmean');
        %         curl_cs = cat(1, curl_cs, zeros(1,size(curl_cs,2)));
        %         curl_cs = cat(2, curl_cs, zeros(size(curl_cs,1),1));
        %         curl_cs(1,:) = 0; curl_cs(:,1) = 0;
        
        %---------------------------------------------------
        % Helmholtz decomp of Flux on fine grid
        % NO need to change grid cell factors for different fluxes
        %---------------------------------------------------
        fprintf('\nHelmholtz decomp...\n\n');
        % BC in fuD&fuR are in 0's
        [phi, psi, fuD, fvD, fuR, fvR] = helmholtz_decomp(fu, fv,...
            onesH, dxCu, dxCv, dyCu, dyCv, [0 0]);
        
        %-- reconstruct DIV & ROT from phi and psi (BCs are all 0's)
        % del_phi = dphi/dx
        [fuD_re, fvD_re] = calc_GxGy_CG(phi,onesH,dxCu,dyCu,dxCv,dyCv,0);
        fuD_re([1 end],:) = 0; fvD_re(:,[1 end]) = 0;
        % z X del_psi = (-dpsi/dy, dpsi/dx)
        [fuR_re, fvR_re] = psi2rot(psi, dyCu, dxCv);
        %
        fprintf( 'RMS(fuD_re-fuD) / RMS(fuD) = %.1e / %.1e \n', rms(fuD_re(:)-fuD(:),'omitnan'), rms(fuD(:),'omitnan') );
        fprintf( 'RMS(fvD_re-fvD) / RMS(fvD) = %.1e / %.1e \n', rms(fvD_re(:)-fvD(:),'omitnan'), rms(fvD(:),'omitnan') );
        fprintf( 'RMS(fuR_re-fuR) / RMS(fuR) = %.1e / %.1e \n', rms(fuR_re(:)-fuR(:),'omitnan'), rms(fuR(:),'omitnan') );
        fprintf( 'RMS(fvR_re-fvR) / RMS(fvR) = %.1e / %.1e \n', rms(fvR_re(:)-fvR(:),'omitnan'), rms(fvR(:),'omitnan') );
        
        %-- Check the accuracy of Helmholtz decomp
        % div & curl
        % [grid.dyCu,grid.dxCv] or [onesU,onesV]
        div_D = calc_div_CG(fuD,fvD,dyCu,dxCv,dxT,dyT,1);
        curl_R = calc_curl_CG(fuR,fvR,dxBu,dyBu);
        curl_D = calc_curl_CG(fuD,fvD,dxBu,dyBu);
        div_R = calc_div_CG(fuR,fvR,dyCu,dxCv,dxT,dyT,1);
        % print the stats
        DF_curl = curl_R - curl;
        DF_div = div_D - divF;
        fprintf( 'RMS (curl_DIV, div_ROT) = %.1e, %.1e \n', rms(curl_D(:),'omitnan'),rms(div_R(:),'omitnan')  );
        fprintf( 'RMS(curl-curl_ROT) / RMS(curl) = %.1e / %.1e \n', rms(DF_curl(:),'omitnan'), rms(curl(:),'omitnan') );
        fprintf( 'RMS(div-div_DIV) / RMS(div) = %.1e / %.1e \n', rms(DF_div(:),'omitnan'), rms(divF(:),'omitnan') );
        
        %---------------------------------------------------
        %       CS-grain psi  --> ROT comp for Flx_CS
        %---------------------------------------------------
        fprintf('\nCS-graining PSI...\n\n');
        % [cs p-grid]
%         psi_cs = sepblockfun(psi(1:end-1,1:end-1),[win_len,win_len],'nanmean');
        psi_cs = imresize(psi,[grid_cs.nih,grid_cs.njh],'Method','nearest'); 

        % back onto cs q-grid, by adding BCs, which must ensure the zero
        % norm flux for the ROT comp from psi_cs
        psi_cs = cat(1, psi_cs, zeros(1,size(psi_cs,2)));
        psi_cs = cat(2, psi_cs, zeros(size(psi_cs,1),1));
        psi_cs(1,:) = 0; psi_cs(:,1) = 0;
        
        %-- check: -z X del{psi_cs} must be non-divergent!
        [fuR_psi_cs, fvR_psi_cs] = psi2rot(psi_cs, grid_cs.dyCu, grid_cs.dxCv);
        div_Rcs = calc_div_CG(fuR_psi_cs,fvR_psi_cs,grid_cs.dyCu,grid_cs.dxCv,grid_cs.dxT,grid_cs.dyT,1);
        fprintf('RMS div{curl{psi_cs}} = %.1e \n',rms(div_Rcs(:)) );
        
        %----------------------------------------------------------------
        %  Get the potential (phi_cs) of <div{flx}> or -<dh/dt>
        %  I.e., solving Poisson eqn
        %                  del2{phi_cs} = SRC = <div{flx}> or -<dh/dt>
        %    with zero norm-flux BC.
        % Note this implies an assumption that SRC can be written in a
        %  form of flux divergence.
        %----------------------------------------------------------------
        fprintf('\nSolving for the DIV comp of <divF>...\n\n');
        src = divF; % div [1/s]    divF = calc_div_CG(fu_do,fv_do,grid.dyCu,grid.dxCv,grid.dxT,grid.dyT,1);
%         src_cs = sepblockfun(src,[win_len,win_len],'nanmean');
        src_cs = imresize(src,[grid_cs.nih,grid_cs.njh],'Method','nearest'); 
        fprintf( 'sum(src_cs)/rms(src_cs) = %.1e / %.1e \n', sum(src_cs(:)), rms(src_cs(:)) );
        
        % PHI ~ [src_cs*A/L], L ~ [dyCu * h / dxCu]
        [~,fuD_phi_cs, fvD_phi_cs] = solve_Poisson(src_cs, areaTcs, onesHcs,...
            grid_cs.dxCu,grid_cs.dyCu,grid_cs.dxCv,grid_cs.dyCv,0);
        
        %-- check if div{fuD_phi_cs} == src_cs
        div_fuD_phi_cs = calc_div_CG(fuD_phi_cs,fvD_phi_cs,grid_cs.dyCu,grid_cs.dxCv,grid_cs.dxT,grid_cs.dyT,1);
        fprintf('RMS(div_fuDcs-src_cs) = %.1e \n',rms(div_fuD_phi_cs- src_cs,'all') );
        %-- check if del_phi is purely divergent
        curl_delphi = calc_curl_CG(fuD_phi_cs,fvD_phi_cs,grid_cs.dxBu,grid_cs.dyBu);
        fprintf('RMS curl_del{phi} = %.1e \n',rms(curl_delphi(:),'omitnan'));
        
        %-----------------------------------------------------------------
        %   CS-grained flux is F_cs = del{phi_cs} - z X del{psi_cs}
        %------------------------------------------------------------------
        % uh
        fu_cs = fuD_phi_cs + fuR_psi_cs;
        fv_cs = fvD_phi_cs + fvR_psi_cs;
        
        %-- check the divergence [fu/m ~ m/s]
        div_F_cs = calc_div_CG(fu_cs,fv_cs,grid_cs.dyCu,grid_cs.dxCv,grid_cs.dxT,grid_cs.dyT,1);
        %     curl_F_cs = calc_curl_CG(fu_cs,fv_cs,grid_cs.dxBu,grid_cs.dyBu);
        
        fprintf( 'RMS(div_F_cs-src_cs) / RMS(src_cs) = %.1e / %.1e \n', ...
            rms(div_F_cs(:)-src_cs(:)), rms(src_cs(:),'omitnan') );
        
    else % not decomp
        
        fu_cs = imresize(fu,[grid_cs.niu,grid_cs.nju],'Method','nearest');
        fv_cs = imresize(fv,[grid_cs.niv,grid_cs.njv],'Method','nearest');
        % enforce no-norm flux BC for u/v
        fu_cs([1 end],:) = 0;
        fv_cs(:,[1 end]) = 0;
        
    end
    
    %-----------------------------------------------------------
    %    assign (multiply the grid len) uh*L ~ m3/s
    %-----------------------------------------------------------
    fu3d_cs(:,:,ik) = fu_cs .* grid_cs.dyCu;
    fv3d_cs(:,:,ik) = fv_cs .* grid_cs.dxCv;
    
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

%% 
function [u_rot, v_rot] = psi2rot(psi, dyCu, dxCv)
% psi on q-; dyCu on u- and dxCv on v-. 
% z X del{psi} = (u_rot, vrot) = (-dpsi_dy, dpsi_dx)
% 
u_rot = - (psi(:,2:end) - psi(:,1:end-1)) ./ dyCu;
v_rot = (psi(2:end,:) - psi(1:end-1,:)) ./ dxCv;
end

end
