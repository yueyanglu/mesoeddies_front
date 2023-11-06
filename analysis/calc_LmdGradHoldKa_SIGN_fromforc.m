% 
% Fit the the along-grad component of Chi, only use one tracer!
% 
% EF = ka (known) - lambda * |del<c>| * SIGN(dC+/dy)
% where C+ is the tracer solution in the "no-EF" runs
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

ntr = 1; 
carry = 9; % 11 for cfc11, 1~8 for idl, 9 for PT

solORref = 1; % 1 for sol, 2 for ref

%-------- dir for  
tmStr = '_180d';
efStr = '_prog_addreld10'; % '_prog', '_prog_addreld10', 
eforc_varnm = 'cforc_rel'; % cforc_rel, cforc_m2,...

% 
grid_dir = [workdir '/MOM6_exp/mom_lowRES/'];  
exp_dir = [workdir '/MOM6_OUT/tr_off_64' tmStr '/'];
h_dir = [workdir '/MOM6_OUT/forc_uvh_64/sol_h' tmStr '/'];

eforc_dir = [exp_dir  'cforc' efStr '/'];
eforc_pt_dir = [exp_dir  'cforc_pt' efStr '/']; 
eforc_cfc_dir = [exp_dir  'cforc_cfcrel125d' efStr '/'];

tr_varnm = 'tr';
% 
if solORref == 1
    disp('Using tr sol to calc params!');
    usetrStr = '_usesol';
    tr_dir = [exp_dir 'sols_idl/tr_forcprog_addrel10d_norel/'];
    tr_pt_dir = [exp_dir 'sols_pt/tr_forcprog_addrel10d_ptrel/'];
    tr_cfc_dir = [exp_dir 'sols_cfc11/tr_forcprog_addrel10d_norel/'];
elseif solORref == 2
    disp('Using ref tr (from HR) to calc params!');
    usetrStr = '_useref';
    tr_dir = [workdir '/MOM6_OUT/tr_off_64/tr_CS/'];
    tr_pt_dir = [workdir '/MOM6_OUT/tr_off_64/tr_pt_CS/']; 
    tr_cfc_dir = [workdir '/MOM6_OUT/tr_off_64/tr_cfc_rel125d_CS/']; 
end
%
tr_noEF_dir = [exp_dir 'sols_idl/tr_noforc/'];
tr_noEF_pt_dir = [exp_dir 'sols_pt/tr_noforc_ptrel/'];
tr_noEF_cfc_dir = [exp_dir 'sols_cfc11/tr_noforc/'];
save_root_dir = [exp_dir 'params/'];

% ------ time
yr_s = 21;
[day_s, day_e, dt] = deal(1, 365*4-6/24, 6/24);  %6/24

t_al = day_s:dt:day_e;
nt_al = length(t_al);
d1yr = 365;
fprintf(1,'Will do time = Y%s D%s.\n',num2str(yr_s), mat2str(t_al));

%-------- read model grid
[grid, ~, ~] = read_grid_MOM([grid_dir '']); % SOLUTION/
nk = length(grid.Layer);
gridsz = size(grid.dxT,1);

%--------------------------------------- how to calc trac grad
cflg = 0; % del_c := h*del_c/dx [c]
if cflg == 0
    fprintf(1,'C-grad NOT multiplied by cell len [c]\n');
elseif cflg == 1
    fprintf(1,'C-grad multiplied by cell len [c*m]\n');
end

%--------------------------------------- hold ka as a known var
k11 = 400;
k22 = 400;
kaStr = '400';

%% --------------------------------------- loop over time

% parfor (icomb = 1, M)%:ncomb 
% for icomb = 1
    
fprintf(1,'Using tracer: %s ...\n',mat2str(carry));

%------------------------dir for saving ------------------------
save_dir = [save_root_dir 'lmddotn_SIGNxmean_ka' kaStr usetrStr efStr ...
     '/C' num2str(carry,'%02d') '/'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
fprintf(1,'Lambda will be saved to: %s\n',save_dir);
    
M = 12;

% ------------------------ calc Chi ------------------------
parfor (it = 1:nt_al, M)
% for it = 1:nt_al

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
    tr_pt_fnm = [tr_pt_dir 'tr__' yrstr '_' dystr '_' hrstr '.nc'];
    eforc_pt_fnm = [eforc_pt_dir 'forc__' yrstr '_' dystr '_' hrstr '.nc'];
    tr_cfc_fnm = [tr_cfc_dir 'tr__' yrstr '_' dystr '_' hrstr '.nc'];
    eforc_cfc_fnm = [eforc_cfc_dir 'forc__' yrstr '_' dystr '_' hrstr '.nc'];
    %
    tr_noEF_fnm = [tr_noEF_dir 'tr__' yrstr '_' dystr '_' hrstr '.nc'];
    tr_noEF_pt_fnm = [tr_noEF_pt_dir 'tr__' yrstr '_' dystr '_' hrstr '.nc'];
    tr_noEF_cfc_fnm = [tr_noEF_cfc_dir 'tr__' yrstr '_' dystr '_' hrstr '.nc'];
    %
    ds_h = ncstruct(h_fnm);
    fprintf(1,'\nReading <h> from: %s...\n', h_fnm);
    %
    if carry == 9
        ds_tr_do = ncstruct(tr_pt_fnm);
        ds_eforc_do = ncstruct(eforc_pt_fnm);
        ds_tr_noEF = ncstruct(tr_noEF_pt_fnm);
        fprintf(1,'Read <c_PT> from: %s\n', tr_pt_fnm);
        fprintf(1,'Read <c_PT> no EF from: %s\n', tr_noEF_pt_fnm);
        fprintf(1,'Read EF for c_PT from: %s\n', eforc_pt_fnm);
    elseif carry>=1 && carry<=8
        ds_tr_do = ncstruct(tr_fnm);
        ds_eforc_do = ncstruct(eforc_fnm);
        ds_tr_noEF = ncstruct(tr_noEF_fnm);
        fprintf(1,'Read <c_IDL> from: %s\n', tr_fnm);
        fprintf(1,'Read <c_IDL> no EF from: %s\n', tr_noEF_fnm);
        fprintf(1,'Read eddy forc from: %s\n', eforc_fnm);
    elseif carry == 11
        ds_tr_do = ncstruct(tr_cfc_fnm);
        ds_eforc_do = ncstruct(eforc_cfc_fnm);
        ds_tr_noEF = ncstruct(tr_noEF_cfc_fnm);
        fprintf(1,'Read <cfc-11> from: %s\n', tr_cfc_fnm);
        fprintf(1,'Read <cfc-11> no EF from: %s\n', tr_noEF_cfc_fnm);
        fprintf(1,'Read eddy forc from: %s\n', eforc_cfc_fnm);
    end

    %----------------------- each layer
    [lmddotn_z] = deal(zeros(grid.nih,grid.njh,nk));

    for ik = 1:nk
        fprintf(1,'\nDoing Z%02d...\n', ik);
        h = ds_h.h(:,:,ik);

        %---------- del<c> and EF
%         [hcx_p,hcy_p, hcxx,hcyy] = deal(zeros(grid.nih,grid.njh));
%         eforc = zeros(grid.nih,grid.njh);
        wichtr = carry;
        if wichtr == 9 || wichtr == 11
            wichtr = 1;
        end
        tr_name = [tr_varnm num2str(wichtr)];
        ef_name = [eforc_varnm num2str(wichtr)];
        fprintf(1,'Read %s.\n', tr_name);
        fprintf(1,'Read %s.\n', ef_name);

        %---------- read in EF and c
        % EF
        eforc = ds_eforc_do.(ef_name)(:,:,ik); % 2d
        % <c>
        tr = ds_tr_do.(tr_name)(:,:,ik);
        % <c> from no-EF run
        tr_noEF = ds_tr_noEF.(tr_name)(:,:,ik);

        %----------  del<c> in [c] & del2<c> in [c/m]
        % delC := h*dc/dx on u/v [c]
        [hcx,hcy] = calc_GxGy_CG(tr,h,grid.dxCu,grid.dyCu,grid.dxCv,grid.dyCv,cflg);
        % BC to 0
        hcx([1 end],:) = 0; hcy(:,[1 end]) = 0;
        % p-grid
        [hcx_p,hcy_p] = uv2p_CG(hcx,hcy);
        % norm of tracer grad
        cgrad_norm = sqrt(hcx_p.^2 + hcy_p.^2);
        % 2nd derivatives: d(hcx)/dx [c/m]
        hcxx = (hcx(2:end,:) - hcx(1:end-1,:)) ./ grid.dxT;
        hcyy = (hcy(:,2:end) - hcy(:,1:end-1)) ./ grid.dyT;

        %--------- no EF tracer grad
        [hcx_noEF,hcy_noEF] = calc_GxGy_CG(tr_noEF,h,grid.dxCu,grid.dyCu,grid.dxCv,grid.dyCv,cflg);
        % BC to 0
        hcx_noEF([1 end],:) = 0; hcy_noEF(:,[1 end]) = 0;
        % p-grid
        [hcx_noEF_p,hcy_noEF_p] = uv2p_CG(hcx_noEF,hcy_noEF);

        %--------- SIGN
        hcy_xmean = mean(hcy_p,1);
        sign2d = sign( repmat(hcy_xmean, [grid.nih 1]) );

        %----------------------- calc lambda [m/s]
        % EF = K11 * cxx + K22 * cyy  - lambda * del<c>
        eforc_ka = k11 .* hcxx + k22 .* hcyy;
        lmddotn = (eforc_ka - eforc) ./ cgrad_norm .* sign2d;
        %
        fprintf('\nCalc Lambda w/ SIGN at Z%02d...\n', ik);
        % lmddotn = filter_extreme(lmddotn,1,99);
        lmddotn_z(:,:,ik) = lmddotn;
    end % k

    %------------------------ save nc
    dim_name = {'xh','yh','zl','Time'};
    dim_length = [grid.nih, grid.njh, nk, 1];
    varname = {'chidotn'};
    dimNum_of_var = {[1,2,3,4]};
    data = {lmddotn_z};

    global_att  = ['eforc = -Chi*|delC|*sign(cy) [m/s]; ' ...
        'carry='  num2str(carry,'%02d') '; ' ...
        'tr_fnm=' tr_fnm, '; h_fnm=' h_fnm '; ' ...
        'tr_pt_fnm=' tr_pt_fnm, '; eforc_pt_fnm=' eforc_pt_fnm '; ' ...
        'eforc_fnm=' eforc_fnm '; eforc_varnm=' eforc_varnm];
    FUN_nc_easywrite_enhanced( savename, dim_name, dim_length,...
        varname, dimNum_of_var, data, global_att )
    fprintf(1,'lambda saved to: %s\n\n',savename);
end % it

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