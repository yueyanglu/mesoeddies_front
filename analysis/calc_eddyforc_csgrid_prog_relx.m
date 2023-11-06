function calc_eddyforc_csgrid_prog_PT_relx(lp)
% 
% Calc eddy forcing: EF = d(<ch> - <c><h>)/dt + div{<U><c>}-<divUc> ... 
%   on CS grid from the pre-saved/-calculated terms.
% 
% The eddy forcing shows up in the evolution eqn of <c>:
% d(<c><h>)dt + div{<U><c>} = EF
% 
% Needs: 
% 1. c, h, dc/dt, dh/dt on HR ---> <c> (ref), d<ch>/dt on LR grid
% 2. U_L on LR
% 3. h_L and dh_L/dt on LR 
% 

homedir = getenv('HOME');
workdir = getenv('WORK');
scradir = getenv('SCRATCH');

addpath(genpath([homedir '/work_MOM']));
addpath(genpath([homedir '/mytoolbox']));
addpath(genpath([homedir '/MyFuncs']));

%% 
% lp = 1;

%-------------------------------------- times
yr_s = 24;
% 
[day_s, day_e, dt_save] = deal(1, 731, 6/24);  
t_al_ful = day_s:dt_save:day_e;
nt_al_ful = length(t_al_ful);
% 
[i_s, dt, i_e] = deal(1, 1600, nt_al_ful);
loops = cell( 1, round((i_e - i_s) / dt) + 1 );
for ii = 1:size(loops,2)
    [is, ie] = deal(i_s + dt * (ii-1), i_s + dt * ii - 1);
    if ie > nt_al_ful
        ie = nt_al_ful;
    end
    loops{ii} = t_al_ful(is:ie);
end
clearvars i_s dt i_e ii is ie
% 
if lp == 0
    t_al = t_al_ful;
else
    t_al = loops{lp};
end
nt_al = length(t_al);
fprintf(1,'This segment do times (nt_al=%d) from D%3.2f to D%3.2f in YR%d.\n',...
    nt_al, t_al(1), t_al(end), yr_s);

% other params
cs_len = 16;
khtr = 100;
method_flx = '4th_order';

%-------------------------------------- dirs
exp_dir = [scradir '/mom_ptemp/'];  
uflx_str = 'prog';
tmStr = '_180d';

uv_dir = [exp_dir 'sol_' uflx_str '/'];  
uv_file_preStr = 'prog__';  nameu = 'uh'; namev = 'vh';
UL_dir = [workdir '/MOM6_OUT/forc_uvh_64/' uflx_str '_CS_decomp' tmStr '/'];
% 
h_dir = [exp_dir 'sol_h/']; % h dir sol_h
hL_dir = [workdir '/MOM6_OUT/forc_uvh_64/sol_h' tmStr '/']; % for prog, use CS HR h (not sol) or sol from uvhm
h_file_preStr = 'h_snap__'; 
% 
% ----- which tracer
carry_al = 1; 
tr_dir = [scradir '/mom_ptemp/sol_tr_pt/'];
save_dir = [workdir '/MOM6_OUT/tr_off_64' tmStr '/cforc_pt_' uflx_str '/' ];

%----- ATM relaxation rate for PT [1/s]
r_relx_atm = 2.9e-8;  
% filename for the relaxation profile 
rel_HR_fnm = [workdir '/MOM6_OUT/forc_uvh_sm17/ICs/trac_init_pt.nc'];
rel_cs_fnm = [workdir '/MOM6_OUT/forc_uvh_64/ICs/trac_init_pt_cs.nc'];
c_ref_HR = ncread(rel_HR_fnm,'tr1');
c_ref_cs = ncread(rel_cs_fnm,'tr1');

% ----
ntr = numel(carry_al);
fprintf(1,'Calc EF for tracer: %s...\n',mat2str(carry_al));
tr_file_preStr = 'tr__';
save_file_preStr = 'forc__'; 

if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%-------------------------------------- grids
%-------- fine model grid
[grid, ~, ~] = read_grid_MOM([exp_dir '']); % SOLUTION/
nk = length(grid.Layer);

%-------- coarse grid, based on "win_len"
grid_cs = build_grid_MOM(grid.nih/cs_len,grid.njh/cs_len,grid.lonq([1 end]),grid.latq([1 end]));
fprintf(1,'Fine grid %d*%d cs-grained to %d*%d.\n',grid.nih, grid.njh, grid_cs.nih, grid_cs.njh);

% dxyU = grid_cs.dxCu .* grid_cs.dyCu;
% dxyV = grid_cs.dxCv .* grid_cs.dyCv;
% dxyT = grid_cs.dxT .* grid_cs.dyT;
onesU = ones(grid.niu,grid.nju);
onesV = ones(grid.niv,grid.njv);
onesUcs = ones(grid_cs.niu,grid_cs.nju);
onesVcs = ones(grid_cs.niv,grid_cs.njv);

%%
%
% Eddy forc at t1 (saved nc name) = "tend at t1" + "prog at t1"
% 

% t_al consists of all the t1
M = 12;
parfor (it = 1:nt_al, M)
% for it = 1:nt_al

    %--------- times: [t1, t1.5, t2]
    t1 = t_al(it);
    % for tend
    [yrstr, dystr, hrstr] = get_timestr(t1, yr_s); 
    % saved NC name
    [yrstr_sv, dystr_sv, hrstr_sv] = get_timestr(t1, yr_s); 
    fprintf(1,'\n time = Y%s-D%s-H%s (it=%d) of %d snapshots ...\n',...
        yrstr,dystr,hrstr,it,nt_al);
    
    %----------------------- check save file
    savefnm = [save_dir save_file_preStr yrstr_sv '_' dystr_sv '_' hrstr_sv '.nc'];
    if exist(savefnm,'file')
        disp(['NC file already exists, skip : ', savefnm]);
        continue
    end
    
    %----------------------- read fine-grid uh, h and tracer nc files
    uv_fnm = [uv_dir uv_file_preStr yrstr '_' dystr '_' hrstr '.nc'];
    UL_fnm = [UL_dir uv_file_preStr yrstr '_' dystr '_' hrstr '.nc'];
    if ~exist(uv_fnm,'file') || ~exist(UL_fnm,'file')
        disp(['NC file NOT exists, skip : ', UL_fnm]);
        continue
    end
    h_fnm = [h_dir h_file_preStr yrstr '_' dystr '_' hrstr '.nc'];
    hL_fnm = [hL_dir h_file_preStr yrstr '_' dystr '_' hrstr '.nc'];
    trac_fnm = [tr_dir tr_file_preStr yrstr '_' dystr '_' hrstr '.nc'];
    %
    ds_uv = ncstruct(uv_fnm);
    ds_csflx = ncstruct(UL_fnm);
    ds_h =  ncstruct(h_fnm);
    ds_hL =  ncstruct(hL_fnm);
    ds_tr = ncstruct(trac_fnm);
    %
    fprintf(1,'\nReading uv from: %s...\n', uv_fnm);
    fprintf(1,'\nReading csg uh from: %s...\n', UL_fnm);
    fprintf(1,'\nReading h from: %s...\n', h_fnm);
    fprintf(1,'\nReading h_L from: %s...\n', hL_fnm);
    fprintf(1,'\nReading tracer from: %s...\n', trac_fnm);
    
    %----------------------- loop tr
    [trforc_tend, trforc_div, trforc_dfus, trforc_rel, trforc] = deal(cell(ntr,1));
    [trforc_tend(:), trforc_div(:), trforc_dfus(:), trforc_rel(:), trforc(:)] = deal({zeros(grid_cs.nih,grid_cs.njh,nk)});
    trforc_mean = cell(ntr,1); trforc_mean(:) = {zeros(grid_cs.nih,grid_cs.njh,nk)};

    for itr = 1:ntr
        wichtr = carry_al(itr);
        tr_name = ['tr' num2str(wichtr)];
        trten_name = ['tr' num2str(wichtr) '_tendency'];
        
        for ik = 1:nk
            
            fprintf(1,'\nDoing Z%02d...\n', ik);
            
            %----- HR flds
            c = ds_tr.(tr_name)(:,:,ik);
            dcdt = ds_tr.(trten_name)(:,:,ik);
            h = ds_h.h(:,:,ik);
            dhdt = ds_h.dhdt(:,:,ik);
            dchdt = c.*dhdt + dcdt.*h;
            % uhL
            [u, v] = deal( ds_uv.u(:,:,ik), ds_uv.v(:,:,ik) );
            [uh, vh] = calc_TFluxes_CG(u.*grid.dyCu, v.*grid.dxCv, h, method_flx);
            % tracer flx [m3/s*c]
            [uc, vc] = calc_TFluxes_CG(uh, vh, c, method_flx);
            div_uc = calc_div_CG(uc, vc,onesU,onesV,grid.dxT,grid.dyT,1);
            % tracer grad: L*h*dc/dx [c*m] & diffus flx: k* L*h*dc/dx [m3/s*c]
            [lhcx,lhcy] = calc_GxGy_CG(c,h,grid.dxCu,grid.dyCu,grid.dxCv,grid.dyCv,1);
            [dfu, dfv] = deal(khtr*lhcx, khtr*lhcy);
            % div_diffus [m/s*c]
            div_dfus = calc_div_CG(dfu, dfv,onesU,onesV,grid.dxT,grid.dyT,1);
            % rh(c*-c)
            rforc = r_relx_atm* h.*(c_ref_HR(:,:,ik) - c);

            %----- ref flds
            c_cs = sepblockfun(c,[cs_len,cs_len],'nanmean');
            dcdt_cs = sepblockfun(dcdt,[cs_len,cs_len],'nanmean');
            dchdt_cs = sepblockfun(dchdt,[cs_len,cs_len],'nanmean');
            %
            uh_L = ds_csflx.(nameu)(:,:,ik); % CS-grained uhl
            vh_L = ds_csflx.(namev)(:,:,ik);
            hL = ds_hL.h(:,:,ik);
            dhLdt = ds_hL.dhdt(:,:,ik);
            
            %--- 1. for tendency
            % d(h_L * <c>)/dt
            dcshLdt = c_cs.*dhLdt + dcdt_cs.*hL;

            %--- 2. for div
            % <div{Uc}>
            div_uc_cs = sepblockfun(div_uc,[cs_len,cs_len],'nanmean');
            % U_L*<c> and div{U_L*<c>}
            [uscs, vscs] = calc_TFluxes_CG(uh_L, vh_L, c_cs, method_flx);
            div_uscs = calc_div_CG(uscs, vscs,onesUcs,onesVcs,grid_cs.dxT,grid_cs.dyT,1);
%             % <c>div<U>
%             cdivu = c_cs .* calc_div_CG(uh_cs, vh_cs,onesUcs,onesVcs,grid_cs.dxT,grid_cs.dyT,1);
%             % <U>*del<c> = div{<U><c>} - <c>div<U>
%             udelc = div_uscs - cdivu;

            %--- 3. for diffusive term
            % <div{k*h*delc}>
            div_dfus_cs = sepblockfun(div_dfus,[cs_len,cs_len],'nanmean');
            % div{k*hL*del<c>}
            [lhcx_cs,lhcy_cs] = calc_GxGy_CG(c_cs,hL,grid_cs.dxCu,grid_cs.dyCu,grid_cs.dxCv,grid_cs.dyCv,1);
            [dfu_cs, dfv_cs] = deal(khtr*lhcx_cs, khtr*lhcy_cs);
            div_dfus_cs2 = calc_div_CG(dfu_cs, dfv_cs,onesUcs,onesVcs,grid_cs.dxT,grid_cs.dyT,1);
           
            %--- 4. for relaxation 
            % <rh(c*-c)> 
            rforc_cs = sepblockfun(rforc,[cs_len,cs_len],'nanmean'); 
            % r*hL*(<c*> - <c>)
            rforc_LR = r_relx_atm * hL.*(c_ref_cs(:,:,ik) - c_cs);

            %----- forcing: tendency + div
            % d(<c><h>)/dt - d<ch>/dt
            trforc_tend{itr}(:,:,ik) = dcshLdt - dchdt_cs;
            % div{<U><c>} - <div_Uc>
            trforc_div{itr}(:,:,ik) = div_uscs - div_uc_cs;
            % <div{k*h*delc}> - div{k*<h>*del<c>}
            trforc_dfus{itr}(:,:,ik) = div_dfus_cs - div_dfus_cs2;
            % <rh(c*-c)>  - r*hL*(<c*> - <c>)
            trforc_rel{itr}(:,:,ik) = rforc_cs - rforc_LR;
            %
            trforc{itr}(:,:,ik) = trforc_tend{itr}(:,:,ik) + ...
                trforc_div{itr}(:,:,ik) + trforc_dfus{itr}(:,:,ik) + ...
                trforc_rel{itr}(:,:,ik);
            %----- d(<c><h>)/dt + div{<U><c>} - div{k*<h>*del<c>}
            trforc_mean{itr}(:,:,ik) = dcshLdt + div_uscs - div_dfus_cs2 ...
                - rforc_LR;
        end % ik

    end % itr
    
    %--------- save to nc
    dim_name = {'xh','yh','zl','Time'};
    dim_length = [grid_cs.nih, grid_cs.njh, nk, 1];
    % t1~ntr, div1~ntr, ...
    varname = { cellstr( num2str(carry_al(:), 'cforc%d') ),...
        cellstr( num2str(carry_al(:), 'cforc_t%d') ),...
        cellstr( num2str(carry_al(:), 'cforc_d%d') ),...
        cellstr( num2str(carry_al(:), 'cforc_dfus%d') ),...
        cellstr( num2str(carry_al(:), 'cforc_r%d') ),...
        cellstr( num2str(carry_al(:), 'cforc_m%d') ) };
    varname = cat(1, varname{:});
    % t1~ntr, div1~ntr, ...
    data = {trforc, trforc_tend, trforc_div, trforc_dfus, trforc_rel, trforc_mean};
    data = cat(1, data{:});
    %
    dimNum_of_var = cell(size(data)); 
    dimNum_of_var(:) = {[1,2,3,4]};
    %
    global_att  = [ 'eforc for CS tracers; ' ...
        'uv_fnm=' uv_fnm '; csflx_fnm=' UL_fnm '; ' ...
        'h_fnm=' h_fnm '; hL_fnm=' hL_fnm '; tr_fnm=' trac_fnm '; ' ...
        'rel_HR_fnm=' rel_HR_fnm '; rel_cs_fnm=' rel_cs_fnm ...
        '; r_relx_atm=' num2str(r_relx_atm) '; ' ...
        'khtr=' num2str(khtr) '; cs_len=' num2str(cs_len)];
    FUN_nc_easywrite_enhanced( savefnm, dim_name, dim_length,...
        varname, dimNum_of_var, data, global_att )
    fprintf(1,'\nCS-grained flds saved to: %s...\n\n', savefnm);
end
delete(gcp('nocreate'));


%% plot 
%{
ik = 1; f_do = trforc_tend{1}(:,:,ik);
hd = pcolor(f_do'); set(hd,'EdgeColor','none'); axis square;
caxis([-1 1]*2e-4);
cmocean('balance'); colorbar
%}
