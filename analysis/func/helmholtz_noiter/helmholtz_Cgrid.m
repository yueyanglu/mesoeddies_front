function [chi, psi, psi_full, urot, vrot, udiv, vdiv] = helmholtz_Cgrid(U, V, GOOD, e1u, e1v, e2u, e2v, WRAP, VERBOSE)
% 
%HELMHOLTZ_CGRID  Solve (u,v) = grad[chi] - z X grad[psi] on a curvilinear Arakawa C grid with holes.
%
% INPUT:
%   U [m by n matrix]: x component of vector, specified at (i-1/2,j).
%   V [m by n matrix]: y component of vector, specified at (i,j-1/2).
%   GOOD [m by n logical matrix]: true where (i,j) is inside the domain.
%   e1u [m by n matrix] dxCu: distance between (i,j) and (i-1,j)
%   e1v [m by n matrix] dxCv: distance across south face of cell 
%   e2u [m by n matrix] dyCu: distance across west face of cell 
%   e2v [m by n matrix] dyCv: distance between (i,j) and (i,j-1)
%   WRAP [logical 2 element vector]: whether to wrap in each dimension
%   VERBOSE 
% 
% OUTPUT:
%   chi [m by n matrix]: potential, on the (i,j) grid
%   psi [m by n matrix]: streamfunction, on the (i-1/2, j-1/2) grid
%   urot [m by n matrix]: rotational part of U, on (i-1/2,j) grid. urot(i,j) = (psi(i,j+1) - psi(i,j)) ./ e2u(i,j);
%   vrot [m by n matrix]: rotational part of V, on (i,j-1/2) grid. vrot(i,j) = -(psi(i+1,j) - psi(i,j)) ./ e1v(i,j);
%   udiv [m by n matrix]: divergent  part of U, on (i-1/2,j) grid. udiv(i,j) = +(chi(i,j) - chi(i-1,j)) ./ e1u(i,j);
%   vdiv [m by n matrix]: divergent  part of V, on (i,j-1/2) grid. vdiv(i,j) = +(chi(i,j) - chi(i,j-1)) ./ e2v(i,j);
%
% Geoff Stanley
% g.stanley@unsw.edu.au // geoffstanley@gmail.com
% v1.0.0  09/06/2017  - initial version
% v1.0.1  20/05/2020  - doco and code cleaning
if nargin < 9 || isempty(VERBOSE)
    VERBOSE = false;
end
if WRAP(1)
    im1 = @(F) circshift(F, [+1 0]);
    ip1 = @(F) circshift(F, [-1 0]);
else
    im1 = @(F) subsasgn(circshift(F, [+1 0]), struct('type','()','subs',{{1,':'}}), 0); %#ok<SUBSASGN>
    ip1 = @(F) subsasgn(circshift(F, [-1 0]), struct('type','()','subs',{{size(F,1),':'}}), 0); %#ok<SUBSASGN>
end
if WRAP(2)
    jm1 = @(F) circshift(F, [0 +1]);
    jp1 = @(F) circshift(F, [0 -1]);
else
    jm1 = @(F) subsasgn(circshift(F, [0 +1]), struct('type','()','subs',{{':',1}}), 0); %#ok<SUBSASGN>
    jp1 = @(F) subsasgn(circshift(F, [0 -1]), struct('type','()','subs',{{':',size(F,2)}}), 0); %#ok<SUBSASGN>
end
% --- Add barriers if data is not periodic
[m,n] = size(GOOD);
% Automatic expansion
e1u = repmat(e1u, m/size(e1u,1), n/size(e1u,2)); 
e1v = repmat(e1v, m/size(e1v,1), n/size(e1v,2)); 
e2u = repmat(e2u, m/size(e2u,1), n/size(e2u,2)); 
e2v = repmat(e2v, m/size(e2v,1), n/size(e2v,2)); 
if WRAP(1) && WRAP(2)
    good_chi = GOOD;
    u = U;
    v = V;
elseif WRAP(1) % so WRAP(2) == false
    n = n + 2;
    good_chi = false(m,n); good_chi(:,2:end-1) = GOOD;
    u = zeros(m,n);              u(:,2:end-1) = U;
    v = zeros(m,n);              v(:,3:end-1) = V(:,2:end);
elseif WRAP(2) % so WRAP(1) == false
    m = m + 2;
    good_chi = false(m,n); good_chi(2:end-1,:) = GOOD;
    u = zeros(m,n);              u(3:end-1,:) = U(2:end,:);
    v = zeros(m,n);              v(2:end-1,:) = V;
else % so WRAP(1) == false and WRAP(2) == false
    n = n + 2;
    m = m + 2;
    good_chi = false(m,n); good_chi(2:end-1,2:end-1) = GOOD;
    u = zeros(m,n);              u(3:end-1,2:end-1) = U(2:end,:);
    v = zeros(m,n);              v(2:end-1,3:end-1) = V(:,2:end);
end
if ~WRAP(1)
    e1u = cat(1, e1u(1,:), e1u, e1u(end,:));
    e1v = cat(1, e1v(1,:), e1v, e1v(end,:));
    e2u = cat(1, e2u(1,:), e2u, e2u(end,:));
    e2v = cat(1, e2v(1,:), e2v, e2v(end,:));
end
if ~WRAP(2)
    e1u = cat(2, e1u(:,1), e1u, e1u(:,end));
    e1v = cat(2, e1v(:,1), e1v, e1v(:,end));
    e2u = cat(2, e2u(:,1), e2u, e2u(:,end));
    e2v = cat(2, e2v(:,1), e2v, e2v(:,end));
end
% --- Select just the largest connected region:
CCocean = CC2periodic(bwconncomp(good_chi,4), WRAP, 'CC');
[~,L] = max(cellfun('length', CCocean.PixelIdxList));
good_chi = false(m,n);
good_chi(CCocean.PixelIdxList{L}) = true;
%--- Set up masks and maps to linear indices, for the different variables
ind.chi = reshape(cumsum(good_chi(:)), m, n);
len.chi = ind.chi(end);
ind.chi(~good_chi) = 0;
good_u = good_chi & im1(good_chi);
good_v = good_chi & jm1(good_chi);
good_psi = good_u & jm1(good_u);
ind.u = reshape(cumsum(good_u(:)), m, n);
len.u = ind.u(end);
ind.u(~good_u) = 0;
ind.v = reshape(cumsum(good_v(:)), m, n);
len.v = ind.v(end);
ind.v = ind.v + len.u;
ind.v(~good_v) = 0;
% --- Find all islands:
CC = CC2periodic(bwconncomp(~good_chi, 8), WRAP, 'CC');
% --- Set one pixel in good_psi to true for each island.
ind_psi_coast = cellfun(@(c) c(end), CC.PixelIdxList);
assert(~any(good_psi(ind_psi_coast)), 'Some border psi pixels are already ''good''...?!');
good_psi(ind_psi_coast) = true;
ind.psi = reshape(cumsum(good_psi(:)), m, n);
len.psi = ind.psi(end);
ind.psi = len.chi + ind.psi;
% --- Fill ind.psi to be uniform around each island
ind_all = reshape(1:m*n, m, n);
ind_all_ip1 = ip1(ind_all);
ind_all_jp1 = jp1(ind_all);
ind_all_ip1jp1 = ip1(ind_all_jp1);
for l = 1:CC.NumObjects
    island_chi = CC.PixelIdxList{l};
    island_psi = [ind_all(island_chi), ind_all_ip1(island_chi), ind_all_jp1(island_chi), ind_all_ip1jp1(island_chi)];
    island_psi(island_psi==0) = [];
    ind.psi(island_psi) = ind.psi(ind_psi_coast(l));
end
clear ind_all ind_all_ip1 ind_all_jp1 ind_all_ip1jp1
if VERBOSE
    fprintf('len.u+len.v=%d, len.chi+len.psi=%d\n', len.u+len.v, len.chi+len.psi);
end
if WRAP(1) && WRAP(2)
    assert(len.chi + len.psi == len.u + len.v, ...
        'Mismatched lengths! len.u+len.v=%d, len.chi+len.psi=%d', len.u+len.v, len.chi+len.psi);
else
    assert(len.chi + len.psi == len.u + len.v + 2, ...
        'Mismatched lengths! len.u+len.v=%d, len.chi+len.psi=%d', len.u+len.v, len.chi+len.psi);
end
% --- Show the grid, to see what we are actually doing:
%
% [xpsi,ypsi] = ndgrid( (1:ind.psi(1))-.5, (1:ind.psi(2))-.5 );
% [xchi,ychi] = ndgrid( (1:ind.chi(1)), (1:ind.chi(2)) );
% [xu,yu] = ndgrid((1:ind.u(1))-.5, (1:ind.u(2)));
% [xv,yv] = ndgrid((1:ind.v(1)), (1:ind.v(2))-.5);
% figure; hold on;
% plot(xchi(good_chi), ychi(good_chi), 'ok');
% plot(xchi(~good_chi), ychi(~good_chi), 's', 'Color', [1 1 1]*.75);
% plot(xpsi(good_psi), ypsi(good_psi), 'xb');
% plot(xpsi(ind_psi_coast), ypsi(ind_psi_coast), 'xg');
% plot(xu(good_u), yu(good_u), '>', 'MarkerSize', 3, 'Color', [1 .5 .5]);
% plot(xv(good_v), yv(good_v), '^', 'MarkerSize', 3, 'Color', [.5 1 .5]);
%
% --- Build matrix:
r = cell(8,1);
c = cell(8,1);
a = cell(8,1);
% Equations for u:
r{1} = ind.u(good_u);
c{1} = ind.chi(good_u);
a{1} = 1./e1u(good_u);
r{2} = r{1};
c{2} = subsref(im1(ind.chi), struct('type', '()', 'subs', {{good_u}}));
a{2} = -a{1};
r{3} = r{1};
c{3} = subsref(jp1(ind.psi), struct('type', '()', 'subs', {{good_u}}));
a{3} = 1./e2u(good_u);
r{4} = r{1};
c{4} = ind.psi(good_u);
a{4} = -a{3};
% Equations for v:
r{5} = ind.v(good_v);
c{5} = ind.chi(good_v);
a{5} = 1./e2v(good_v);
r{6} = r{5};
c{6} = subsref(jm1(ind.chi), struct('type', '()', 'subs', {{good_v}}));
a{6} = -a{5};
r{7} = r{5};
c{7} = ind.psi(good_v);
a{7} = 1./e1v(good_v);
r{8} = r{5};
c{8} = subsref(ip1(ind.psi), struct('type', '()', 'subs', {{good_v}}));
a{8} = -a{7};
szA = len.chi + len.psi;
if WRAP(1) && WRAP(2)
    % We have: len.u + len.v == len.chi + len.psi == size(A,1) == size(A,2)
    % but rank(full(A)) == size(A,1) - 2, so A is rank deficient by 2.
    % Remove one equation for chi and one for psi,
    % and add one equation for chi and one for psi that specify their values at a particular spot.
    r = cell2mat(r);
    c = cell2mat(c);
    a = cell2mat(a);
    A = sparse(r, c, a, szA, szA);
    A(1,:) = 0; A(1,1) = 1;
    A(len.chi+1,:) = 0; A(len.chi+1,len.chi+1) = 1;
    y = [u(good_u); v(good_v);];
    y([1 len.chi+1]) = 0;
    
else
    % We have: len.u + len.v + 2 == len.chi + len.psi == size(A,1) == size(A,2)
    % but rank(full(A)) == size(A,1) - 2, so A is rank deficient by 2.
    % Add one equation for chi and one for psi that specify their values at a particular spot.
    r = [cell2mat(r); szA-1; szA];
    c = [cell2mat(c); 1; len.chi + 1];
    a = [cell2mat(a); 1; 1];
    A = sparse(r, c, a, szA, szA);
    y = [u(good_u); v(good_v); 0; 0];
end
% --- Solve linear system (using mldivide):
mytic = tic; 
x = A \ y;
if VERBOSE
    fprintf('Solved by mldivide in time %.4f, rms error = %.4g\n', toc(mytic), rms(y - A*x));
end
% --- Construct chi and psi:
chi = zeros(m,n);
chi(good_chi) = x(1:len.chi);
psi = x(ind.psi);
% --- Construct rotational and divergent velocities
if nargout >= 3
    urot = (jp1(psi) - psi) ./ e2u;
    urot = trim(urot);
end
if nargout >= 4
    vrot = -(ip1(psi) - psi) ./ e1v;
    vrot = trim(vrot);
end
if nargout >= 5
    udiv = (chi - im1(chi)) ./ e1u;
    udiv(~good_chi | im1(~good_chi)) = 0;
    udiv = trim(udiv);
end
if nargout >= 6
    vdiv = (chi - jm1(chi)) ./ e2v;
    vdiv(~good_chi | jm1(~good_chi)) = 0;
    vdiv = trim(vdiv);
end

% only works when WRAP = [0 0]!!!
psi_full = psi(2:end,2:end);

% original
chi = trim(chi);
psi = trim(psi);

    function var = trim(var)
        if ~WRAP(1) && ~WRAP(2)
            var = var(2:end-1,2:end-1);
        elseif ~WRAP(1)
            var = var(2:end-1,:);
        elseif ~WRAP(2)
            var = var(:,2:end-1);
            % else, leave var alone
        end
    end
end