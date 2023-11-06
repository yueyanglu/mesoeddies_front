
function [grid, stat, fname] = read_grid_MOM(dir)
% 
% Read the geometry, vertical coordinate, depth and ocean status nc files
% 
% e.g. dir = '/glade/scratch/yueyanglu/momexp_fine_60s/SOLUTION/'
% Syntax: [grid, stat] = read_grid_MOM(dir)
% 
%----
geom_fnm = [dir 'ocean_geometry.nc'];
vert_fnm = [dir 'Vertical_coordinate.nc'];
dep_fnm  = [dir 'Depth_list.nc'];
stat_fnm = [dir 'ocean.stats.nc'];

%---- read as struct
geom = ncstruct(geom_fnm);
vert = ncstruct(vert_fnm);
dep  = ncstruct(dep_fnm);
stat = ncstruct(stat_fnm);

%---- combine
grid = catstruct(geom, vert, dep);

%---- add some params that are not provided by MOM6
Lx = grid.lonq(end) - grid.lonq(1);   % [km]
Ly = grid.latq(end) - grid.latq(1); 
dx = grid.dxT(1); %  consisten with dxT
dy = grid.dyT(2); 
[nih, njh] = size(grid.geolon);
[niu, nju] = deal(nih+1, njh);
[niv, njv] = deal(nih, njh+1);
[niq, njq] = deal(nih+1, njh+1);

%---- save to grid
grid.Lx = Lx; grid.Ly = Ly;  
grid.dx = dx; grid.dy = dy; 
grid.nih = nih; grid.njh = njh; 
grid.niu = niu; grid.nju = nju; 
grid.niv = niv; grid.njv = njv; 
grid.niq = niq; grid.njq = njq; 

% ---- save file names
fname.dep_fnm = dep_fnm;
fname.geom_fnm = geom_fnm;
fname.stat_fnm = stat_fnm;
fname.vert_fnm = vert_fnm;



