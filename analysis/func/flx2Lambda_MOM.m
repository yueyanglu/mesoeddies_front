function [lmdu,lmdv] = flx2Lambda_MOM(lhs,cxp,cyp)
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
%      lmd  [nih,njh]: [m/s]
% 
% lhs = lambda * del<c>
% 

%% check size of input

% must be 3D
if numel(size(lhs)) >  3
    error('f2K function accepts 3D fields only (j-i-ndist) !!');
end

% inputs must have the same size
if ~isequal(numel(lhs),numel(cxp),numel(cyp))
    error('Size of the inputs must be the same !!!');
end

% size
[ni,nj,ntr] = size(lhs);

ntr_min = 2;
if ntr < ntr_min
    error(['Must have at least ' num2str(ntr_min) ' traceres !!!']);
end

%% calc K

[lmdu, lmdv] = deal(NaN * zeros(ni,nj)); 

% Note this loop will miss the N and E boundaries of the NE
% corner bin.
% ~1min for one snapshot (sequel)
disp('doing L only')
for i = 1:ni
    for j = 1:nj
        
        %------------------------------------ LHS, column vec [div1;div2;div3;..]
        F = squeeze(lhs(i,j,:)); 
        
        %------------------------------------ matrix of del2c & grad
        % dCdx & dCdy on p [ntr-1]
        Gxp = squeeze(cxp(i,j,:));
        Gyp = squeeze(cyp(i,j,:));
        % Form the LHS matrix
        G = [Gxp, Gyp];
        
        %------------------------------------ inverse, -G * KL = D
        if all(~isnan(G),'all') && all(~isnan(F),'all')
            
            if rank(G) == ntr_min
%                 tol = norm(G)/1e2;
%                 KL = lsqminnorm(G,F);
%                 KL = G \ F;
                KL = pinv(G) * F;
%                 KL = G \ F;  % K = - G \ F; or K =  pinv(G) * F
            else
                KL(1:ntr_min) = NaN;
            end
            
            lmdu(i,j) = KL(1); lmdv(i,j) = KL(2);
        end
    end
end
