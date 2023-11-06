function [fu_sm,fv_sm] = smooth_uv(grid,win_len,fu,fv,ifdecomp,ifverb)
%[fu_sm,fv_sm] = smooth_uv(grid,win_len,fu,fv,ifverb)
% This function smooth a 2D vector field (e.g., uh and vh)!
% 
% Input: 
%    fu2d - u*h [m2/s]
%    ifdecomp - "0" for simple filtering; "1" for decomp-based
% 
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
if nargin <= 5
    ifverb = 0;
    disp(['ifverb  = ', num2str(ifverb)])
end
if ifdecomp == 0
    ifverb = 0;
end

onesH = ones(grid.nih,grid.njh);
[dxCu, dxCv, dyCu, dyCv] = deal(grid.dxCu, grid.dxCv, grid.dyCu, grid.dyCv);
[dxBu, dyBu, dxT, dyT] = deal(grid.dxBu, grid.dyBu, grid.dxT, grid.dyT);
dxyT = grid.dxT .* grid.dyT;
dxyU = grid.dxCu .* grid.dyCu;
dxyV = grid.dxCv .* grid.dyCv;
% dxyB = grid.dxBu .* grid.dyBu;
%
ifpad = 1;
% [xTrunc, yTrunc] = deal( (win_len+1)/2 );

%% do the snapshot


if ifdecomp 
    
    %---------------------------------------------------
    % decompose-based filtering
    %---------------------------------------------------
    % div [fu/m] and curl [fu/m] of flux
    divF = calc_div_CG(fu,fv,dyCu,dxCv,dxT,dyT,1);
    curl = calc_curl_CG(fu,fv,dxBu,dyBu);
    
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
    if ifverb
        fprintf( 'RMS(fuD_re-fuD) / RMS(fuD) = %.1e / %.1e \n', rms(fuD_re(:)-fuD(:),'omitnan'), rms(fuD(:),'omitnan') );
        fprintf( 'RMS(fvD_re-fvD) / RMS(fvD) = %.1e / %.1e \n', rms(fvD_re(:)-fvD(:),'omitnan'), rms(fvD(:),'omitnan') );
        fprintf( 'RMS(fuR_re-fuR) / RMS(fuR) = %.1e / %.1e \n', rms(fuR_re(:)-fuR(:),'omitnan'), rms(fuR(:),'omitnan') );
        fprintf( 'RMS(fvR_re-fvR) / RMS(fvR) = %.1e / %.1e \n', rms(fvR_re(:)-fvR(:),'omitnan'), rms(fvR(:),'omitnan') );
    end
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
    if ifverb
        fprintf( 'RMS (curl_DIV, div_ROT) = %.1e, %.1e \n', rms(curl_D(:),'omitnan'),rms(div_R(:),'omitnan')  );
        fprintf( 'RMS(curl-curl_ROT) / RMS(curl) = %.1e / %.1e \n', rms(DF_curl(:),'omitnan'), rms(curl(:),'omitnan') );
        fprintf( 'RMS(div-div_DIV) / RMS(div) = %.1e / %.1e \n', rms(DF_div(:),'omitnan'), rms(divF(:),'omitnan') );
    end
    %---------------------------------------------------
    %       SMOOTH psi  --> ROT comp for Flx_SM
    %---------------------------------------------------
    fprintf('\nSMOOTHING PSI...\n\n');
    % [ p-grid]
    %     psi_sm = smooth_geom_CG(psi,dxyB,win_len,win_len,2,ifpad);
    %     psi_sm([1 end],:) = 0; psi_sm(:,[1 end]) = 0;
    
    % back onto cs q-grid, by adding BCs, which must ensure the zero
    % norm flux for the ROT comp from psi_cs
    psi_sm = smooth_geom_CG(psi(1:end-1,1:end-1),dxyT,win_len,win_len,2,ifpad);
    psi_sm = cat(1, psi_sm, zeros(1,size(psi_sm,2)));
    psi_sm = cat(2, psi_sm, zeros(size(psi_sm,1),1));
    psi_sm(1,:) = 0; psi_sm(:,1) = 0;
    
    %-- check: -z X del{psi_cs} must be non-divergent!
    [fuR_psi_sm, fvR_psi_sm] = psi2rot(psi_sm, grid.dyCu, grid.dxCv);
    div_Rcs = calc_div_CG(fuR_psi_sm,fvR_psi_sm,grid.dyCu,grid.dxCv,grid.dxT,grid.dyT,1);
    if ifverb
        fprintf('RMS div{curl{psi_cs}} = %.1e \n',rms(div_Rcs(:)) );
    end
    %----------------------------------------------------------------
    %  Get the potential (phi_cs) of <div{flx}> or -<dh/dt>
    %  I.e., solving Poisson eqn
    %                  del2{phi_cs} = SRC = <div{flx}> or -<dh/dt>
    %    with zero norm-flux BC.
    % Note this implies an assumption that SRC can be written in a
    %  form of flux divergence.
    %----------------------------------------------------------------
    fprintf('\nSolving for the DIV comp of <divF>...\n\n');
    src = divF; % div [1/s]
    src_sm = smooth_geom_CG(src,dxyT,win_len,win_len,2,ifpad);
    
    % PHI ~ [src_cs*A/L], L ~ [dyCu * h / dxCu]
    % del_PHI ~ [PHI/dx]
    [~,fuD_phi_sm, fvD_phi_sm] = solve_Poisson(src_sm, dxyT, onesH,...
        grid.dxCu,grid.dyCu,grid.dxCv,grid.dyCv,0,0);
    
    %-- check if div{fuD_phi_cs} == src_cs
    div_fuD_phi_sm = calc_div_CG(fuD_phi_sm,fvD_phi_sm,grid.dyCu,grid.dxCv,grid.dxT,grid.dyT,1);
    %-- check if del_phi is purely divergent
    curl_delphi = calc_curl_CG(fuD_phi_sm,fvD_phi_sm,grid.dxBu,grid.dyBu);
    
    if ifverb
        fprintf( 'sum(src_cs)/rms(src_cs) = %.1e / %.1e \n', sum(src_sm(:)), rms(src_sm(:)) );
        fprintf('RMS(div_fuDcs-src_cs) = %.1e \n',rms(div_fuD_phi_sm- src_sm,'all') );
        fprintf('RMS curl_del{phi} = %.1e \n',rms(curl_delphi(:),'omitnan'));
    end
    %-----------------------------------------------------------------
    %   CS-grained flux is F_cs = del{phi_cs} - z X del{psi_cs}
    %------------------------------------------------------------------
    % uh
    fu_sm = fuD_phi_sm + fuR_psi_sm;
    fv_sm = fvD_phi_sm + fvR_psi_sm;
    
    %-----------------------------------------------------------------
    % if ifcat
    %
    %     % for near bdry: use deecompossed
    %     [fu_decomp, fv_decomp] = deal(fu_sm,fv_sm);
    %     fu_decomp(1+xTrunc:end-xTrunc,1+yTrunc:end-yTrunc) = 0;
    %     fv_decomp(1+xTrunc:end-xTrunc,1+yTrunc:end-yTrunc) = 0;
    %
    %     % for inside: use simply filtered
    %     fu_filt = smooth_geom_CG(fu,dxyU,win_len,win_len,2,ifpad);
    %     fv_filt = smooth_geom_CG(fv,dxyV,win_len,win_len,2,ifpad);
    %     fu_filt([1:xTrunc end-xTrunc+1:end],:) = 0; fu_filt(:,[1:yTrunc end-yTrunc+1:end]) = 0;
    %     fv_filt([1:xTrunc end-xTrunc+1:end],:) = 0; fv_filt(:,[1:yTrunc end-yTrunc+1:end]) = 0;
    %
    %     % near the bdry- decomp; inside - just filter
    %     fu_sm = fu_decomp + fu_filt;
    %     fv_sm = fv_decomp + fv_filt;
    % end
    
    
    %-- check the divergence [fu/m ~ m/s]
    div_F_sm = calc_div_CG(fu_sm,fv_sm,grid.dyCu,grid.dxCv,grid.dxT,grid.dyT,1);
    %     curl_F_cs = calc_curl_CG(fu_cs,fv_cs,grid_cs.dxBu,grid_cs.dyBu);
    
    if ifverb
        fprintf( 'RMS(div_F_sm-src_cs) / RMS(src_cs) = %.1e / %.1e \n', ...
            rms(div_F_sm(:)-src_sm(:)), rms(src_sm(:),'omitnan') );
    end
    
else
    
    fu_sm = smooth_geom_CG(fu, dxyU, win_len, win_len, 2, ifpad);
    fv_sm = smooth_geom_CG(fv, dxyV, win_len, win_len, 2, ifpad);
    
end

%% 
function [u_rot, v_rot] = psi2rot(psi, dyCu, dxCv)
% psi on q-; dyCu on u- and dxCv on v-. 
% z X del{psi} = (u_rot, vrot) = (-dpsi_dy, dpsi_dx)
% 
u_rot = - (psi(:,2:end) - psi(:,1:end-1)) ./ dyCu;
v_rot = (psi(2:end,:) - psi(1:end-1,:)) ./ dxCv;
end

end
