function L = coef_Laplacian_CG(h,dxCu,dyCu,dxCv,dyCv)
% 
% The coefficients for the 2D Laplacian/Poission equation, del2_PHI = f,
% where PHI is the potential, using 5-points convention.
% All fields are defined on C-grid (MOM6 sym). PHI and f is assumed to be column
% vectors with a size of nih*njh-by-1 (i first).
% [(1,1), (2,1),..., (ni,1), ... ]'
% This is consistent with 'calc_GxGy_CG'!
%
% L = coef_Laplacian_CG(ones(5,5),ones(6,5),ones(6,5),ones(5,6),ones(5,6));
% 


[nih,njh] = size(h);
[niu,~] = size(dxCu);
[~,njv] = size(dyCv);
if nih ~= niu-1 || njh ~= njv-1 
    error('Fields must be defined on C-grid with [x,y] !!!');
end

onemu = 1.e-12; % thickness of vanishing layer

% layer thickness at the interface by harmonic mean [u-grid]
h_harm_u = harmon( max(h(1:end-1,:),onemu), max(h(2:end,:),onemu) ); 
h_harm_u = cat(1, h_harm_u(1,:), h_harm_u, h_harm_u(end,:)); % cp BC points

% layer thickness at the interface by harmonic mean [v-grid]
h_harm_v = harmon( max(h(:,1:end-1),onemu), max(h(:,2:end),onemu) ); 
h_harm_v = cat(2, h_harm_v(:,1), h_harm_v, h_harm_v(:,end));

% dxCu = ones(nih+1,nih);
% dyCv = ones(nih,nih+1);
% dxT = ones(nih,nih);
% dyT = ones(nih,nih);

%%
% L is a nih*njh-by-nih*njh matrix
[I,J,L1d] = deal([]);

mm = 0; % linear index of non-zeros entries in L

% interior grid points
for j = 1:njh %1:njh
    for i = 1:nih
        
        %--- local coef
        % i and j in dxCu(i,j) should change with the sym/asym C-Grid !!!
        ci_m1 = dyCu(i,j) * h_harm_u(i,j) / dxCu(i,j);   % PHI(i-1,j) coefficient
        ci_p1 = dyCu(i+1,j) * h_harm_u(i+1,j) / dxCu(i+1,j); % PHI(i+1,j) 
        cj_m1 = dxCv(i,j) * h_harm_v(i,j)/ dyCv(i,j);   % PHI(i,j-1)  coefficient
        cj_p1 = dxCv(i,j+1) * h_harm_v(i,j+1)/ dyCv(i,j+1); % PHI(i,j+1) coefficient
        cij = - (ci_m1 + ci_p1 + cj_m1 + cj_p1); % PHI(i,j) coefficient
        
%         ci_m1 = 1 / dxCu(i,j) / dxT(i,j);   % PHI(i-1,j) coefficient
%         ci_p1 = 1 / dxCu(i+1,j) / dxT(i,j); % PHI(i+1,j) coefficient
%         cj_m1 = 1 / dyCv(i,j) / dyT(i,j);   % PHI(i,j-1)  coefficient
%         cj_p1 = 1 / dyCv(i,j+1) / dyT(i,j); % PHI(i,j+1) coefficient
%         cij = - (ci_m1 + ci_p1 + cj_m1 + cj_p1); % PHI(i,j) coefficient
        
        %--- linear index of PHI col vec, also the col index for L
        colL = sub2ind([nih,njh], i,j);
        
        %--- PHI(i,j)
        mm = mm + 1;
        mm0 = mm; % position of the point on the main diagonal
        rowL = colL;
        [I(mm), J(mm)] = deal(rowL,colL);
        L1d(mm) = cij;
        
        %--- PHI(i-1,j) WEST
        if i > 1 
            mm = mm + 1;
            rowL = colL-1; % & the row index for L
            [I(mm), J(mm)] = deal(rowL,colL);
            L1d(mm) = ci_m1;
        elseif i == 1 % zero flux BC
            L1d(mm0) = L1d(mm0) + ci_m1;
        end
       
        %--- PHI(i+1,j) EAST
        if i < nih
            mm = mm + 1;
            rowL = colL+1;
            [I(mm), J(mm)] = deal(rowL,colL);
            L1d(mm) = ci_p1;
        elseif i == nih
            L1d(mm0) = L1d(mm0) + ci_p1;
        end
        
        %--- PHI(i,j-1)
        if j > 1
            mm = mm + 1;
            rowL = colL-nih;
            [I(mm), J(mm)] = deal(rowL,colL);
            L1d(mm) = cj_m1;
        elseif j == 1
            L1d(mm0) = L1d(mm0) + cj_m1;
        end
        
        %--- PHI(i,j+1)
        if j < njh
            mm = mm + 1;
            rowL = colL+nih;
            [I(mm), J(mm)] = deal(rowL,colL);
            L1d(mm) = cj_p1;
        elseif j == njh
            L1d(mm0) = L1d(mm0) + cj_p1;
        end
        
%         L(colL-1,colL) = ci_m1;    % PHI(i-1,j) 
%         L(colL+1,colL) = ci_p1;    % PHI(i+1,j) 
%         L(colL,colL)   = cij;   % PHI(i,j) 
%         L(colL-nih,colL) = cj_m1;    % PHI(i,j-1) 
%         L(colL+nih,colL) = cj_p1;    % PHI(i,j+1) 
        
    end
end

L = sparse(I,J,L1d, (nih-0)*(njh-0), (nih-0)*(njh-0));  % full(L)

end

function h=harmon(x,y)
h=2*x.*y./(x+y);
end


