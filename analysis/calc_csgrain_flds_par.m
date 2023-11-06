function calc_csgrain_flds_par(lp, h_or_flx)
%Coarse-grain h or uh fields. 
% Run a time segment in a parfors loop
% h_or_flx: 1 for h; 2 for flux

homedir = getenv('HOME');
workdir = getenv('WORK');
scradir = getenv('SCRATCH');

addpath(genpath([homedir '/work_MOM']));
addpath(genpath([homedir '/mytoolbox']));
addpath(genpath([homedir '/MyFuncs']));


%% 
ifdecomp = 1;
dt_flxmean = 6/24; % time-interval over which flux is averaged [day]

%-------------------------------------- times
yr_s = 20;
% 
[day_s, day_e, dt_save] = deal(216, 366.0, 6/24);  % 1/24
t_al_ful = day_s:dt_save:day_e;
nt_al_ful = length(t_al_ful);
% 
[i_s, dt, i_e] = deal(1, 10, nt_al_ful);
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
% fprintf(1,'This segment do time (nt_al=%d) = Y%d D%s.\n',nt_al, yr_s, mat2str(t_al));
fprintf(1,'This segment do times (nt_al=%d) from D%3.2f to D%3.2f in YR%d.\n',...
    nt_al, t_al(1), t_al(end), yr_s);

% other params
fillvalue = 1.0e+20;
cs_len = 16;   
% 
if ifdecomp
    decomStr = 'decomp';
else
    decomStr = 'nodecomp';
end

%-------------------------------------- dirs
exp_dir = [scradir '/mom_ptemp/'];  
flx_dir = [exp_dir 'sol_uvhm/'];   % uv dir
flx_file_preStr = 'uvh_mean__';  nameu = 'uh'; namev = 'vh';
save_varname_fu = nameu; save_varname_fv = namev;
% 
h_dir = [exp_dir 'sol_h/']; % h dir
h_file_preStr = 'h_snap__'; nameh = 'h';
save_varname_h = nameh;

if h_or_flx == 1 % h only
    save_file_preStr = h_file_preStr;
    save_dir = [workdir '/MOM6_OUT/forc_uvh_64/sol_h_CS/'];
elseif h_or_flx == 2 % flx only
    save_file_preStr = flx_file_preStr;
    save_dir = [workdir '/MOM6_OUT/forc_uvh_64/uvhm_CS_' decomStr '/'];
else
    error('Define ''h_or_flx'' to either 1 or 2 !!!!')
end

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

%% time

M = 12;

parfor (it = 1:nt_al, M)
% parfor it = 1:nt_al

    %----------------------- current time
    if h_or_flx == 1
        nday = t_al(it);
    elseif h_or_flx == 2 
        % The NC name of time-mean flux over [t1~t2] is (t1+t2)/2!
        nday = (t_al(it) + t_al(it)+dt_flxmean) / 2;
    end
    [yrstr, dystr, hrstr] = get_timestr(nday, yr_s); % for tend
    fprintf(1,'\n time = Y%s-D%s-H%s (it=%d) of %d snapshots ...\n',...
        yrstr,dystr,hrstr,it,nt_al);
    
    %----------------------- check save file
    savefnm = [save_dir save_file_preStr yrstr '_' dystr '_' hrstr '.nc'];
    if exist(savefnm,'file') 
        disp(['NC file already exists, skip : ', savefnm]);
        continue
%         system(['rm ' savefnm]); warning(['NC file already exists, delete it : ', savefnm]);
    end    
    
    %----------------------- read fields
    if h_or_flx == 1
        
        %--- read
        h_fnm = [h_dir h_file_preStr yrstr '_' dystr '_' hrstr '.nc'];
        h3d =  ncread(h_fnm, nameh);
        dhdt3d =  ncread(h_fnm, 'dhdt');
        fprintf(1,'\nReading %s from: %s...\n', nameh, h_fnm);
        
        %--- calc
        h3d_cs = zeros(grid_cs.nih, grid_cs.njh, nk);
        dhdt3d_cs = zeros(grid_cs.nih, grid_cs.njh, nk);
        for ik = 1:nk
            h3d_cs(:,:,ik) = sepblockfun(h3d(:,:,ik),[cs_len,cs_len],'nanmean');
            dhdt3d_cs(:,:,ik) = sepblockfun(dhdt3d(:,:,ik),[cs_len,cs_len],'nanmean');
        end
        
        %---- save
        nccreate(savefnm,save_varname_h,'Format','netcdf4','Datatype','double',...
            'Dimensions',{'xh',grid_cs.nih,'yh',grid_cs.njh,'zl',nk,'Time',1},...
            'FillValue',fillvalue);
        nccreate(savefnm,'dhdt','Format','netcdf4','Datatype','double',...
            'Dimensions',{'xh',grid_cs.nih,'yh',grid_cs.njh,'zl',nk,'Time',1},...
            'FillValue',fillvalue);
        ncwrite(savefnm,save_varname_h,h3d_cs);
        ncwrite(savefnm,'dhdt',dhdt3d_cs);
        fprintf(1,'\nCS-grained %s saved to: %s...\n\n', nameh, savefnm);
        
    elseif h_or_flx == 2
        
        %--- read
        flx_fnm = [flx_dir flx_file_preStr yrstr '_' dystr '_' hrstr '.nc'];
        fu3d = ncread(flx_fnm, nameu);
        fv3d = ncread(flx_fnm, namev);
        fprintf(1,'\nReading flx %s & %s from: %s...\n', nameu, namev, flx_fnm);
        
        %--- calc
        [fu3d_cs,fv3d_cs] = csgrain_flds(grid,grid_cs,fu3d,fv3d,ifdecomp);
        
        %---- save
        nccreate(savefnm,save_varname_fu,'Format','netcdf4','Datatype','double',...
            'Dimensions',{'xq',grid_cs.niu,'yh',grid_cs.nju,'zl',nk,'Time',1},...
            'FillValue',fillvalue);
        nccreate(savefnm,save_varname_fv,'Format','netcdf4','Datatype','double',...
            'Dimensions',{'xh',grid_cs.niv,'yq',grid_cs.njv,'zl',nk,'Time',1},...
            'FillValue',fillvalue);
        %
        ncwrite(savefnm,save_varname_fu,fu3d_cs);
        ncwrite(savefnm,save_varname_fv,fv3d_cs);
        fprintf(1,'\nCS-grained %s saved to: %s...\n\n', nameu, savefnm);
    end
   
end
delete(gcp('nocreate'));

%----------------------- save important info to a .mat
savemat = [save_dir 'params.mat'];
if ~exist(savemat,'file')
    save(savemat, '*_dir','*_preStr','cs_len','t_al','grid','grid_cs','fillvalue');
    fprintf(1,'\nImportant configures saved to: %s...\n\n', savemat);
end

%% plot to check
% ik = 1;
% fu = fu3d_cs(:,:,ik); fv = fv3d_cs(:,:,ik); div = calc_div_CG(fu,fv,ones(grid_cs.niu,grid_cs.nju),ones(grid_cs.niv,grid_cs.njv),grid_cs.dxT,grid_cs.dyT,1);
% hd = pcolor(div'); set(hd,'EdgeColor','none'); axis square;
% caxis([-1 1]*2e-4);
% cmocean('balance'); colorbar

