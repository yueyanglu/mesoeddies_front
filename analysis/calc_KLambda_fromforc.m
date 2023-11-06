% 
% Fit the kappa-lambda model from the pre-calc eddy forcing
% 
% lhs = - K11 * cxx - 2*K12 * cxy - K22 * cyy + lambda * del<c>
% lhs = d(c'h')/dt + ... + div{U'c'} + ...
% 
% NEED; 
%   Eddy forcing [m/s*c]
%   layer thickness
%   tracer, from which cxx, cyy and delC are calc
% 

clear
homedir = getenv('HOME');
workdir = getenv('WORK');
scradir = getenv('SCRATCH');
addpath(genpath([homedir '/work_MOM']));
addpath(genpath([homedir '/mytoolbox']));
addpath(genpath([homedir '/MyFuncs']));

%% confg

%---- set from shell
% klay = klaySh; %klaySh;
% icomb = icSh; ntr = ncSh;
% wichSM = smSh;
% ifExtra = extraSh;
% ifiso = isoSh; % 1 iso; 0 aniso
% wichTerm = tmSh; % not used for now

ntr = 8; % icomb = '6' for 12357; '12' for 12457
carry_al = [1 2 3 4 5 6 7 8];

% ntr = 5; % ntr=3: icomb =22 for C159; ntr=5 icomb =7 for C12357; =17 for C12457 
% carry_al = [1 2 3 4 5 6 7 8 9];

ifiso = 1; 
solORref = 2; % 1 for sol, 2 for ref

%-------- dir for  
tmStr = '_180d';
efStr = '_prog'; % _uvhm or '_prog', '_prog_addreld10', '_uvhm_addreld10'
% 
grid_dir = [workdir '/MOM6_exp/mom_lowRES/'];  
exp_dir = [workdir '/MOM6_OUT/tr_off_64' tmStr '/'];
h_dir = [workdir '/MOM6_OUT/forc_uvh_64/sol_h' tmStr '/'];
eforc_dir = [exp_dir  'cforc' efStr '/'];
eforc_pt_dir = [exp_dir  'cforc_pt' efStr '/']; 
eforc_varnm = 'cforc'; % cforc_m1, cforc_m2,...
tr_varnm = 'tr';
% 
if solORref == 1
    disp('Using tr sol to calc params!');
    usetrStr = '_usesol';
    tr_dir = [exp_dir 'sols_idl/tr_forc_addreld10_norel/'];
    tr_pt_dir = [exp_dir 'sols_pt/tr_forc_addreld10_ptrel_norel/'];
elseif solORref == 2
    disp('Using ref tr (from HR) to calc params!');
    usetrStr = '_useref';
    tr_dir = [workdir '/MOM6_OUT/tr_off_64/tr_CS/'];
    tr_pt_dir = [workdir '/MOM6_OUT/tr_off_64/tr_pt_CS/']; 
end
save_root_dir = [exp_dir 'params/'];

% ------ time
yr_s = 21;
% [day_s, day_e, dt] = deal(1.0+6/24, 730.0, 6/24); 
[day_s, day_e, dt] = deal(1, 730, 6/24); 

t_al = day_s:dt:day_e;
nt_al = length(t_al);
d1yr = 365;
fprintf(1,'Will do time = Y%s D%s.\n', num2str(yr_s), mat2str(t_al));

%-------- iso model or aniso model
if ifiso
    isoStr = 'iso';
else
    isoStr = 'ani';
end

%-------- read model grid
[grid, ~, ~] = read_grid_MOM([grid_dir '']); % SOLUTION/
nk = length(grid.Layer);
gridsz = size(grid.dxT,1);
fillvalue = 1.0e+20;

%--------------------------------------- how to calc trac grad
cflg = 0; % del_c := h*del_c/dx [c]
if cflg == 0
    fprintf(1,'C-grad NOT multiplied by cell len [c]\n');
elseif cflg == 1
    fprintf(1,'C-grad multiplied by cell len [c*m]\n');
end

%--------------------------------------- all possible tracer sets
trac_comb = nchoosek(carry_al,ntr);
ncomb = size(trac_comb,1);

%% --------------------------------------- loop over tracer set
M = 12;

parfor (icomb = 1, M)%:ncomb 
    
    if icomb == 0
        carries = carry_al; % use all avail tracers to over-determine
    else
        if icomb > size(trac_comb,1)
            warning('Id of combination exceeds # of avail comb, set to the last one')
            carries = trac_comb(end,:);
        else
            carries = trac_comb(icomb,:);
        end
    end
    fprintf(1,'Using tracers: %s ...\n',mat2str(carries));
    
    %------------------------ dir for saving ------------------------
    %
    save_dir = [save_root_dir 'klmd_' isoStr '' usetrStr efStr ...
         '/C' num2str(carries,'%02d') '/'];
    if ~exist(save_dir,'dir')
        mkdir(save_dir);
    end
    fprintf(1,'K&Lambda will be saved to: %s\n',save_dir);
    
    %------------------------ calc k&L ------------------------
    for it = 1:nt_al
        
        %----------------------- current time
        [yrstr, dystr, hrstr] = get_timestr(t_al(it), yr_s);
        fprintf(1,'\n time = Y%s-D%s-H%s (it=%d) of %d snapshots ...\n',...
            yrstr,dystr,hrstr,it,nt_al);
        
        %----------------------- filename to be saved
        savename = [save_dir 'K__' yrstr '_' dystr '_' hrstr '.nc'];
        % if exist, skip iteration
        if exist(savename,'file')
            fprintf(1,'\nK&L already exists, SKIP: %s\n\n',savename);
            continue
        end
        
        %----------------------- read structures of <h>, <c> and EF [m/s*c]
        h_fnm = [h_dir 'h_snap__' yrstr '_' dystr '_' hrstr '.nc'];
        tr_fnm = [tr_dir 'tr__' yrstr '_' dystr '_' hrstr '.nc'];
        eforc_fnm = [eforc_dir 'forc__' yrstr '_' dystr '_' hrstr '.nc'];
        ds_h = ncstruct(h_fnm);
        ds_tr = ncstruct(tr_fnm);
        ds_eforc = ncstruct(eforc_fnm);
        fprintf(1,'\nReading <h> from: %s...\n', h_fnm);
        fprintf(1,'Read <c> from: %s\n', tr_fnm);
        fprintf(1,'Read eddy forc from: %s\n', eforc_fnm);
        % 
        tr_pt_fnm = [tr_pt_dir 'tr__' yrstr '_' dystr '_' hrstr '.nc'];
        eforc_pt_fnm = [eforc_pt_dir 'forc__' yrstr '_' dystr '_' hrstr '.nc'];
        if any(carries == 9)
            ds_tr_pt = ncstruct(tr_pt_fnm);
            ds_eforc_pt = ncstruct(eforc_pt_fnm);
            fprintf(1,'Read <c_PT> from: %s\n', tr_pt_fnm);
            fprintf(1,'Read EF for c_PT from: %s\n', eforc_pt_fnm);
        end

        %----------------------- each layer
        [K11_z,K12_z,K22_z,lmdu_z,lmdv_z] = deal(zeros(grid.nih,grid.njh,nk));
        
        for ik = 1:nk
            fprintf(1,'\nDoing Z%02d...\n', ik);
            h = ds_h.h(:,:,ik);
            
            %---------- put together all del<c> and eforc
            [hcx_p,hcy_p,del2c,hcxx,hcxy,hcyy] = deal(zeros(grid.nih,grid.njh,ntr));
            eforc = zeros(grid.nih,grid.njh,ntr);
            
            for itr = 1:ntr
                wichtr = carries(itr);
                if wichtr == 9
                    wichtr = 1;
                    ds_eforc_do = ds_eforc_pt;
                    ds_tr_do = ds_tr_pt;
                    fprintf(1,'Read PT (tr9).\n');
                else
                    ds_eforc_do = ds_eforc;
                    ds_tr_do = ds_tr;
                end
                tr_name = [tr_varnm num2str(wichtr)];
                ef_name = [eforc_varnm num2str(wichtr)];
                fprintf(1,'Read %s.\n', tr_name);
                fprintf(1,'Read %s.\n', ef_name);
                
                %---------- read in EF and c
                % EF
                eforc(:,:,itr) = ds_eforc_do.(ef_name)(:,:,ik); % 2d
                % <c>
                tr = ds_tr_do.(tr_name)(:,:,ik);
            
                %----------  del<c> in [c] & del2<c> in [c/m]
                % delC := h*dc/dx on u/v [c]
                [hcx,hcy] = calc_GxGy_CG(tr,h,grid.dxCu,grid.dyCu,grid.dxCv,grid.dyCv,cflg);
                % BC to 0
                hcx([1 end],:) = 0; hcy(:,[1 end]) = 0;
                % p-grid
                [hcx_p(:,:,itr),hcy_p(:,:,itr)] = uv2p_CG(hcx,hcy);
                
                % Laplacian: { d(hcx)*L + d(hcy)*L} / dxdy2 [c/m]
                del2c(:,:,itr) = calc_div_CG(hcx,hcy,grid.dyCu,grid.dxCv,grid.dxT,grid.dyT,1);
                % 2nd derivatives: d(hcx)/dx [c/m]
                hcxx(:,:,itr) = (hcx(2:end,:) - hcx(1:end-1,:)) ./ grid.dxT;
                hcyy(:,:,itr) = (hcy(:,2:end) - hcy(:,1:end-1)) ./ grid.dyT;
                % cat cxy along y-direction
                hcx_temp = cat(2, zeros(grid.niu,1), hcx, zeros(grid.niu,1));
                hcxy_q = (hcx_temp(:,2:end) - hcx_temp(:,1:end-1)) ./ grid.dyBu;
                hcxy(:,:,itr) = (hcxy_q(1:end-1,1:end-1) + hcxy_q(1:end-1,2:end) + ...
                    hcxy_q(2:end,2:end) + hcxy_q(2:end,1:end-1)) / 4;
            end % ntr
        
            %----------------------- calc k [m2/s] & lambda [m/s]
            % -eforc = - K11 * cxx - 2*K12 * cxy - K22 * cyy + lambda * del<c>
            lhs = - eforc;
            fprintf('\nCalc K & Lambda at Z%02d...\n', ik);
            [K11,K12,K22,lmdu,lmdv] = flx2KLambda_MOM(lhs,hcxx,hcxy,hcyy,hcx_p,hcy_p,ifiso);
            %
%             K11 = filter_extreme(K11,1,99);
%             K12 = filter_extreme(K12,1,99);
%             K22 = filter_extreme(K22,1,99);
%             lmdu = filter_extreme(lmdu,1,99);
%             lmdv = filter_extreme(lmdv,1,99);
            %
%             K11(isnan(K11)) = 0; lmdu(isnan(lmdu)) = 0; lmdv(isnan(lmdv)) = 0; 
            K11_z(:,:,ik) = K11; K12_z(:,:,ik) = K12; K22_z(:,:,ik) = K22;
            lmdu_z(:,:,ik) = lmdu;  
            lmdv_z(:,:,ik) = lmdv;
        end % k
        
        %------------------------ save nc
        if ifiso
            dim_name = {'xh','yh','zl','Time'};
            dim_length = [grid.nih, grid.njh, nk, 1];
            varname = {'K11', 'lmdu', 'lmdv'};
            dimNum_of_var = {[1,2,3,4], [1,2,3,4], [1,2,3,4]};
            data = {K11_z, lmdu_z, lmdv_z};
        else
            dim_name = {'xh','yh','zl','Time'};
            dim_length = [grid.nih, grid.njh, nk, 1];
            varname = {'K11', 'K12', 'K22', 'lmdu', 'lmdv'};
            dimNum_of_var = {[1,2,3,4], [1,2,3,4], [1,2,3,4]};
            data = {K11_z, K12_z, K22_z, lmdu_z, lmdv_z};
        end
        global_att  = ['-eforc = - k [m2/s] + Chi [m/s] with ' isoStr '; ' ...
            'carries='  num2str(carries,'%02d') '; ' ...
            'tr_fnm=' tr_fnm, '; h_fnm=' h_fnm '; ' ...
            'tr_pt_fnm=' tr_pt_fnm, '; eforc_pt_fnm=' eforc_pt_fnm '; ' ...
            'eforc_fnm=' eforc_fnm '; eforc_varnm=' eforc_varnm];
        FUN_nc_easywrite_enhanced( savename, dim_name, dim_length,...
            varname, dimNum_of_var, data, global_att )
        fprintf(1,'K&lambda saved to: %s\n\n',savename);
        
    end % it

end % icomb
delete(gcp('nocreate'));

%% reconstructed eddy forcing := -kij*del2c + lmd * delc
%{

ntr = length(carries);

[term1, term2] = deal(zeros(grid.nih,grid.njh,ntr));
for itr = 1:ntr
    if ifiso
        term1(:,:,itr) = - K11 .* (hcxx(:,:,itr) + hcyy(:,:,itr));
        term2(:,:,itr) = lmdu .* hcx_p(:,:,itr) + lmdv .* hcy_p(:,:,itr);
    else
        term1(:,:,itr) = - (K11.*hcxx(:,:,itr) + 2*K11.*hcyy(:,:,itr) + K22.*hcyy(:,:,itr))
        term2(:,:,itr) = lmdu .* hcx_p(:,:,itr) + lmdv .* hcy_p(:,:,itr);
    end
end
lhs_re = term1 + term2;

% ---- plot
itr = 3;
% clim = [-1e-4, 1e-4];
clim = [-5e-5, 5e-5];
f_do = {term1(:,:,itr), term2(:,:,itr)}; titles = {'-kdel2c', 'lmd*delc'};
% f_do = {lhs(:,:,itr), lhs_re(:,:,itr)};

figure
for icel = 1:2
    subplot(1,2,icel)
    hd = pcolor(grid.lonh,grid.lath,f_do{icel}');  set(hd, 'EdgeColor', 'none');
    axis square
    cmocean('balance');
    caxis(clim)
    title(titles{icel})
    colorbar
end

%% plot lhs and mean fields used to calc params

itr = 3;
win_extra = 31;
dd = 0;
clim = [-1e-4, 1e-4];
figure
f_do = lhs(:,:,itr);
% f_do = smooth_geom_CG(f_do,dxyT,win_extra,win_extra);
hd = pcolor(f_do(dd+1:end-dd,dd+1:end-dd)');  set(hd, 'EdgeColor', 'none');
axis square
cmocean('balance');
caxis(clim)
title('eddy forcing')
colorbar

% ------
clim = [-1e-8, 1e-8];
figure
f_do = hcyy(:,:,itr);
hd = pcolor(grid.lonh,grid.lath,f_do');  set(hd, 'EdgeColor', 'none');
axis square
cmocean('balance');
caxis(clim)
colorbar

% ------
clim = [-1e-3, 1e-3];
figure
subplot(121)
f_do = hcx_p(:,:,itr);
hd = pcolor(grid.lonh,grid.lath,f_do');  set(hd, 'EdgeColor', 'none');
axis square
cmocean('balance');
caxis(clim)
title('d<c>/dx')
colorbar
% 
subplot(122)
f_do = hcy_p(:,:,itr);
hd = pcolor(grid.lonh,grid.lath,f_do');  set(hd, 'EdgeColor', 'none');
axis square
cmocean('balance');
caxis(clim)
title('d<c>/dy')
colorbar

%% --------- 2D maps of kappa & Lambda
% [grid, ~, grid_files] = read_grid_MOM('/Users/yueyang/Documents/MATLAB/MOM/momexp_fine/'); % SOLUTION/

clim = [-5e3, 5e3; -.2, .2];
% clim = [-2e3, 2e3; -.05, .05].*1e2;

cmname = 'balance';

font = 'DejaVu Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);

if ifiso
    titls = {'\kappa', '\chi_{u}', '\chi_{v}'}; % /\langle{h}\rangle
    plt_fields = {K11, lmdu, lmdv}; % .*dpm
    
    [ha, ~] = tight_subplot(1,3,[.02 .02],[.02 .02],[.06 .05]);
    for icel = 1:3
        axes(ha(icel))
        f_do = plt_fields{icel};
        
        hd = pcolor(f_do');  set(hd, 'EdgeColor', 'none');
        axis square
        cmocean(cmname);
        title([ '$' titls{icel} '$'],'interpreter','latex','FontSize',16);
        set(gca,'xticklabel','','yticklabel','')
        if icel == 1
            
            caxis(clim(1,:))
            cb = colorbar;
            cb.Orientation = 'horizontal';
            posPlt = ha(icel).Position;
            posCb = [posPlt(1)  0.20  posPlt(3)  0.02];
            cb.Position = posCb;
            cb.Title.String = 'm^2{s^{-1}}';
            cb.Title.Position = [80 -30 0];
        else
            caxis(clim(2,:))
            if icel == 2
                cb = colorbar;
                cb.Orientation = 'horizontal';
                posPlt = ha(icel).Position;
                posCb = [posPlt(1)+posPlt(3)/2  0.20  posPlt(3)  0.02];
                cb.Position = posCb;
                cb.Title.String = 'm {s^{-1}}';
                cb.Title.Position = [80 -30 0];
            end
        end
    end
else % not iso
    
    titls = {'\kappa_{11}', '\kappa_{12}', '\kappa_{22}'}; % /\langle{h}\rangle
    plt_fields = {K11, K12, K22}; % .*dpm
    [ha, ~] = tight_subplot(1,3,[.02 .02],[.02 .02],[.06 .05]);
    for icel = 1:3
        axes(ha(icel))
        %     subplot(2,2,icel)
        f_do = plt_fields{icel};
        f_do = smooth_geom_HYCOM(f_do,grid.dxT.*grid.dyT,win_extra,win_extra);
        
        hd = pcolor(grid.lonh,grid.lath,f_do');  set(hd, 'EdgeColor', 'none');
        axis square
        cmocean(cmname);
        title([ '$' titls{icel} '$'],'interpreter','latex','FontSize',16);
        set(gca,'xticklabel','','yticklabel','')
        caxis(clim(1,:))
        cb = colorbar;
        cb.Orientation = 'horizontal';
        posPlt = ha(icel).Position;
        posCb = [posPlt(1)  0.20  posPlt(3)  0.02];
        cb.Position = posCb;
        cb.Title.String = 'm^2{s^{-1}}';
        cb.Title.Position = [80 -30 0];
    end
    
    font = 'DejaVu Sans';
    figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
        titls = {'\chi_{u}', '\chi_{v}'}; % /\langle{h}\rangle
    plt_fields = {lmdu, lmdv}; % .*dpm
    [ha, ~] = tight_subplot(1,2,[.02 .02],[.02 .02],[.06 .05]);
    for icel = 1:2
        axes(ha(icel))
        %     subplot(2,2,icel)
        f_do = plt_fields{icel};
%         f_do = smooth_geom_HYCOM(f_do,grid.dxT.*grid.dyT,win_extra,win_extra);
        
        hd = pcolor(grid.lonh,grid.lath,f_do');  set(hd, 'EdgeColor', 'none');
        axis square
        cmocean(cmname);
        title([ '$' titls{icel} '$'],'interpreter','latex','FontSize',16);
        set(gca,'xticklabel','','yticklabel','')
        caxis(clim(2,:))
        cb = colorbar;
        cb.Orientation = 'horizontal';
        posPlt = ha(icel).Position;
        posCb = [posPlt(1)  0.17  posPlt(3)  0.02];
        cb.Position = posCb;
        cb.Title.String = 'm^2{s^{-1}}';
        cb.Title.Position = [80 -30 0];
    end
end

%% kij and lambda



%%
%
% figure;hd = pcolor(grid.lonh,grid.lath,eddysrc.tend{1}(:,:,itr)'); 
% set(hd, 'EdgeColor', 'none'); caxis([-1e-3 1e-3]); cmocean('balance'); colorbar
% 
% 
figure
for icel = 1:ncel
    subplot(2,2,icel)
    hd = pcolor(grid.lonh,grid.lath,eddysrc.tend{icel}(:,:,itr)');
    set(hd, 'EdgeColor', 'none'); caxis([-5e-4 5e-4]); cmocean('balance'); colorbar
end

figure
for icel = 1:ncel
    subplot(2,2,icel)
    hd = pcolor(grid.lonh,grid.lath,eddysrc.divuc{icel}(:,:,itr)');
    set(hd, 'EdgeColor', 'none'); caxis([-5e-4 5e-4]); cmocean('balance'); colorbar
end

%%
c = sc.cS + sc.cE;
dcdt = sc.dcSdt + sc.dcEdt;
h = sc.hS + sc.hE;
dhdt = sc.dhSdt + sc.dhEdt;
dchdt = c .* dhdt + h .* dcdt;
% 
uc = sc.uscs + sc.usce + sc.uecs + sc.uece;
vc = sc.vscs + sc.vsce + sc.vecs + sc.vece;
div_uc = calc_div_CG(uc, vc, onesU,onesV,grid.dxT,grid.dyT,1); 
% diffusive flux := khtr*  (L *h*delta_c /dx) [c* m3/s]
% del<c> := L*h*dc/dx on u/v [c*m]
[hlcx,hlcy] = calc_GxGy_CG(c,h,grid.dxCu,grid.dyCu,grid.dxCv,grid.dyCv,1);
khtr = 50;
% [c* m3/s]
[dfu, dfv] = deal(khtr*hlcx, khtr*hlcy);
% [c* m/s]
div_diffus = calc_div_CG(dfu, dfv, onesU,onesV,grid.dxT,grid.dyT,1); 

figure
subplot(131)
hd = pcolor(grid.lonh,grid.lath,dchdt');
axis square
set(hd, 'EdgeColor', 'none'); caxis([-5e-4 5e-4]); cmocean('balance'); colorbar
%     
subplot(132)
hd = pcolor(grid.lonh,grid.lath,div_uc');
axis square
set(hd, 'EdgeColor', 'none'); caxis([-5e-4 5e-4]); cmocean('balance'); colorbar
% 
subplot(133)
hd = pcolor(grid.lonh,grid.lath, dchdt' + div_uc' - div_diffus');
axis square
set(hd, 'EdgeColor', 'none'); caxis([-5e-4 5e-4]); cmocean('balance'); colorbar
% delta = (dchdt + div_uc - div_diffus) ./ div_uc;
% mean(abs(delta(:)))
%}