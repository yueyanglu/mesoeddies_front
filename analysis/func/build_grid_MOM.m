function grid = build_grid_MOM(nih,njh,xlim,ylim)
% 
% Build the Cartesian grid. Fields are defined as [nx-ny] on C-grid
% 
% xlim:    Zonal bounds  [km]
% ylim:    Meridional bounds [km]
% 
% Syntax: grid = build_grid_MOM(nih,njh,xlim,ylim)
% 

% -------- params
[niu, nju] = deal(nih+1, njh);
[niv, njv] = deal(nih, njh+1);
[niq, njq] = deal(nih+1, njh+1);
Lx = diff(xlim); % [km] 
Ly = diff(ylim); % [km]

% -------- Coordinates (consistent with the symmetric mode of MOM6) [km]
[lonh,lonq,dx] = FDGrid(xlim(1), xlim(2), nih);    
[lath,latq,dy] = FDGrid(ylim(1), ylim(2), njh);
lonh = lonh'; lonq = lonq'; lath = lath'; latq = latq';
[geolon,geolat] = ndgrid(lonh,lath); % p     
[geolonu,geolatu] = ndgrid(lonq,lath); % u             
[geolonv,geolatv] = ndgrid(lonh,latq); % v         
[geolonb,geolatb] = ndgrid(lonq,latq); % q           

% -------- grid spacing [m]
factor = 1e3;
dx = dx *factor;
dy = dy *factor;
% p-cell
dxT = diff(geolonu,1,1) *factor; 
dyT = diff(geolatv,1,2) *factor;
% u-cell
dxCu = diff(geolon,1,1) *factor; dxCu = cat(1, dxCu(1,:), dxCu, dxCu(end,:)); 
dyCu = diff(geolatb,1,2) *factor; 
% v-cell
dxCv = diff(geolonb,1,1) *factor; 
dyCv = diff(geolat,1,2) *factor; dyCv = cat(2, dyCv(:,1), dyCv, dyCv(:,end));
% q-cell
dxBu = diff(geolonv,1,1) *factor; dxBu = cat(1, dxBu(1,:), dxBu, dxBu(end,:)); 
dyBu = diff(geolatu,1,2) *factor; dyBu = cat(2, dyBu(:,1), dyBu, dyBu(:,end));

% -------- areas for cells [m2]
Ah = dxT.*dyT;
Aq = dxBu.*dyBu;

% -------- save
grid = struct('nih',nih,'njh',njh, 'niu',niu,'nju',nju, 'niv',niv,'njv',njv, ...
    'niq',niq,'njq',njq, 'Lx',Lx,'Ly',Ly, 'dx',dx,'dy',dy, ...
    'lonh',lonh,'lath',lath,'lonq',lonq,'latq',latq,...
    'geolon',geolon,'geolat',geolat,'geolonb',geolonb,'geolatb',geolatb, ...
    'geolonu',geolonu,'geolatu',geolatu,'geolonv',geolonv,'geolatv',geolatv, ...
    'dxT',dxT,'dyT',dyT, 'dxCu',dxCu,'dyCu',dyCu, 'dxCv',dxCv,'dyCv',dyCv, ...
    'dxBu',dxBu,'dyBu',dyBu, 'Ah',Ah, 'Aq',Aq);



