function [tend, divuc, divdf, dfflx, termsStr] = get_eddyforc(means,eddies,ds_h,grid,khtr)
%Calc mean and eddy tracer forcing from the <c>, U'<c> terms, on the same
% grid.
%
% INPUT: 
%   means  - struct containing <c>, d<c>/dt, <U><c>...
%   eddies - struct containing c', dc'/dt, <U>c'...
%   ds_h   - struct containing <h>, h', d<h>/dt, dh'/dt
%   grid   - struct 
%   khtr   - isopycnal tracer diffusivity [m2/s]
% 
% OUTPUT:
%   tend  - cell array containing mean and eddy tracer tendency
%   divuc - tracer flux div
%   divdf - tracer diffusive flux div
% 

termsStr = {'<U><c>', ' U''<c> ', ' <U>c'' ', ' U''c'' ', 'sum'};
ncel = numel(termsStr);

% for convenience
[cS, cE] = deal(means.cs, eddies.ce);
[dcSdt, dcEdt] = deal(means.dcsdt, eddies.dcedt);
[hS, hE, dhSdt, dhEdt] = deal(ds_h.hS, ds_h.hE, ds_h.dhSdt, ds_h.dhEdt);

%--------------------------------------------------
%           tracer (eddy) forcing [m/s*c]
%--------------------------------------------------
[tend, divuc, divdf] = deal( cell(1,ncel) );
[uedflx, vedflx, hcx, hcy, dfu, dfv] = deal( cell(1,ncel) );

%------------ tendency [m/s*c] ----------------
% d(<h><c>)dt, d(h'<c>)dt, d(<h>c')dt, d(h'c')dt, sum_eddy
tend{1} = cS .* dhSdt + dcSdt .* hS;
tend{2} = cS .* dhEdt + dcSdt .* hE;
tend{3} = cE .* dhSdt + dcEdt .* hS;
tend{4} = cE .* dhEdt + dcEdt .* hE;
tend{5} = tend{2} + tend{3} + tend{4};

%------------ eddy flux div [m/s*c] ----------------
% eddy flx [m2/s*c]:  <U><c> U'<c>, <U>c' & U'c'
uedflx{1} = means.uscs;   vedflx{1} = means.vscs;
uedflx{2} = eddies.uecs;  vedflx{2} = eddies.vecs;
uedflx{3} = eddies.usce;  vedflx{3} = eddies.vsce;
uedflx{4} = eddies.uece;  vedflx{4} = eddies.vece;
uedflx{5} = uedflx{2} + uedflx{3} + uedflx{4};
vedflx{5} = vedflx{2} + vedflx{3} + vedflx{4};
%--- div[m/s*c] above, note the eddy flx in uhc [m2/s*c]
for icel = 1:ncel
    divuc{icel} = calc_div_CG(uedflx{icel}, vedflx{icel}, ...
        grid.dyCu,grid.dxCv,grid.dxT,grid.dyT,1);
end

%------------ diffusive flux div [m/s*c] ----------------
% del_c: h*dc/dx [c*m] & diffus flx: k* h*dc/dx [m2/s*c]
[hcx{1},hcy{1}] = calc_GxGy_CG(cS,hS,grid.dxCu,grid.dyCu,grid.dxCv,grid.dyCv,0);
[hcx{2},hcy{2}] = calc_GxGy_CG(cS,hE,grid.dxCu,grid.dyCu,grid.dxCv,grid.dyCv,0);
[hcx{3},hcy{3}] = calc_GxGy_CG(cE,hS,grid.dxCu,grid.dyCu,grid.dxCv,grid.dyCv,0);
[hcx{4},hcy{4}] = calc_GxGy_CG(cE,hE,grid.dxCu,grid.dyCu,grid.dxCv,grid.dyCv,0);
% [hdcdx,hdcdy] = calc_GxGy_CG(cS+cE,hS+hE,grid.dxCu,grid.dyCu,grid.dxCv,grid.dyCv,0);
hcx{5} = hcx{2} + hcx{3} + hcx{4};
hcy{5} = hcy{2} + hcy{3} + hcy{4};
%--- div[m/s*c]
for icel = 1:ncel
    [dfu{icel}, dfv{icel}] = deal(khtr*hcx{icel}, khtr*hcy{icel});
    divdf{icel} = calc_div_CG(dfu{icel}, dfv{icel}, ...
        grid.dyCu,grid.dxCv,grid.dxT,grid.dyT,1);
end

%------------ save diffusive fluxes ----------------
dfflx.dfu = dfu; dfflx.dfv = dfv; dfflx.khtr = khtr;



