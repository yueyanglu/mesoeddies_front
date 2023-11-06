function [K11,K12,K22,lmdu,lmdv] = flx2KLambda_MOM(lhs,cxx,cxy,cyy,cxp,cyp,ifiso)
% 
% 'flx2KLambda_MOM.m' calculates the instantaneous eddy diffusivity 
%  and 'advective vector' from the divergence/advection of eddy tracer flux.
% 
%   Input:
%      lhs   [nih,njh,ntr]: local eddy forcing [c*m/s]      
%      cxx   [nih,njh,ntr]: 2nd derivative in [c/m], e.g. d(hdc/dx)/dx
%      cxy, cyy 
%      cxp   [nih,njh,ntr]: x-grad of c on p-point [c], i.e., h*dc/dx
%      cyp   [nih,njh,ntr]: y-grad of c on p-point [c]
%      ifiso              : if use isotropic diffusion (0 for anisotropic)
%
%   Output:
%      K    [nih,njh]: [m2/s]
%      lmd  [nih,njh]: [m/s]
% 
% lhs = - K11 * cxx - 2*K12 * cxy - K22 * cyy + lambda * del<c>
% 

%% check size of input

% must be 3D
if numel(size(lhs)) >  3
    error('f2K function accepts 3D fields only (j-i-ndist) !!');
end

% inputs must have the same size
if ~isequal(numel(lhs),numel(cxx),numel(cyy),numel(cxp),numel(cyp))
    error('Size of the inputs must be the same !!!');
end

% size
[ni,nj,ntr] = size(lhs);

% iso diffusion: K*laplacian, 3 unknowns, >= 3 tracers;
% anis diffusion: s11*cxx + s12*cxy + s22*cyy, 5 unknowns, >= 5 tracers;
if ifiso
    ntr_min = 3;
else
    ntr_min = 5;
end

if ntr < ntr_min
    error(['Must have at least ' num2str(ntr_min) ' traceres because'...
        ' ifiso = ' num2str(ifiso) ' !!!']);
end

%% calc K

[K11, K12, K22, lmdu, lmdv] = deal(NaN * zeros(ni,nj)); 

% Note this loop will miss the N and E boundaries of the NE
% corner bin.
% ~1min for one snapshot (sequel)
for i = 1:ni
    for j = 1:nj
        
        %------------------------------------ LHS, column vec [div1;div2;div3;..]
        F = squeeze(lhs(i,j,:)); 
        
        %------------------------------------ matrix of del2c & grad
        % cxx, cxy, cyy [ntr-1]
        Gxx = squeeze(cxx(i,j,:));
        Gxy = squeeze(cxy(i,j,:));
        Gyy = squeeze(cyy(i,j,:));
        % dCdx & dCdy on p
        Gxp = squeeze(cxp(i,j,:));
        Gyp = squeeze(cyp(i,j,:));
        
        % Form the LHS matrix
        % if iso, ndist-by-3; if aniso, ndist-by-5
        if ifiso
            G = [-Gxx-Gyy, Gxp, Gyp];
        else
            G = [-Gxx, -2*Gxy, -Gyy, Gxp, Gyp];
        end
        
        %------------------------------------ inverse, -G * KL = D
        if all(~isnan(G),'all') && all(~isnan(F),'all')
            
            %if rank(G) == ntr_min
                %tol = norm(G)/1e2;
                %KL = lsqminnorm(G,F); % tol
		%KL = G \ F;
		KL = pinv(G) * F;  %pinv is equiv to KL = lsqminnorm(G,F);
                %KL = pinv(G) * F;  % K = - G \ F; or K = - pinv(G) * F
            %else
            %    KL(1:ntr_min) = NaN;
            %end
            
            if ifiso
                K11(i,j) = KL(1); K12(i,j) = KL(1); K22(i,j) = KL(1);
                lmdu(i,j) = KL(2); lmdv(i,j) = KL(3);
            else
                K11(i,j) = KL(1); K12(i,j) = KL(2); K22(i,j) = KL(3);
                lmdu(i,j) = KL(4); lmdv(i,j) = KL(5);
            end
        end
    end
end
