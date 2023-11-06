function [Kxx,Kxy,Kyx,Kyy] = flx2K_MOM(dcdx_u,dcdy_v,fu,fv)
% 
% 'flx2K_HYCOM.m' calculates the instantaneous eddy diffusivity tensor 
%  and/or asscociated variables from the given large-scale tracer and eddy
%  tracer flux.
% 
%   Input:
%      cxu  [niu,nju,ntr]: zonal tracer grad on u-point [m*c ~ len*c]
%      cyv  [niv,njv,ntr]: meridional tracer grad on v-point 
%      fu   [niu,nju,ntr]: Zonal eddy tracer flx      [m/s*c ~ len*c]
%      fv   [niv,njv,ntr]: Meridional eddy tracer flx
%
%   Output:
%     Kxx,Kxy   [niu,niu]: m2/s
%     Kyx,Kyy   [niv,niv]: 
% 

%% check size of input

% must be 2D
if numel(size(dcdx_u)) >  3
    error('f2K function accepts 3D fields only (j-i-ndist) !!');
end

% size
[niu,nju,ntr] = size(dcdx_u);
[niv,njv,~] = size(dcdy_v);
[nih,njh] = deal(niv, nju);

if ntr < 2
    error('Must have at least 2 traceres !!!');
end

%% interpolate large-scale tracer gradient on v/u points

dcdx_v = NaN * zeros(niv,njv,ntr);
dcdy_u = NaN * zeros(niu,nju,ntr);

for itr = 1:ntr
    
    % u- and v-grid
    [cxu, cyv] = deal(dcdx_u(:,:,itr), dcdy_v(:,:,itr));
    
    % Find dCdx on the v- points [1:niv,2:njv-1], NS bndry are untouched
    cxv = (cxu(1:end-1,1:end-1) + cxu(1:end-1,2:end) + ...
        cxu(2:end,2:end) + cxu(2:end,1:end-1)) / 4; % from SW, AC
    
    % Find dCdy on the u- points [2:niu-1,1:nju], WE bndry are untouched
    cyu = (cyv(1:end-1,1:end-1) + cyv(1:end-1,2:end)+...
        cyv(2:end,2:end) + cyv(2:end,1:end-1)) / 4;
    
    % ASSIGN
    [dcdx_v(:,2:end-1,itr), dcdy_u(2:end-1,:,itr)]...
        = deal(cxv, cyu);
    
end
  

%% calc K

[Kxx,Kxy] = deal(NaN * zeros(niu,nju)); 
[Kyx,Kyy] = deal(NaN * zeros(niv,njv)); 
% [angu,angv] = deal(NaN * zeros(nj,ni)); 

% 
% Note that MOM6 uses north-east convention centered on the h-points
% This will miss the S and W boundaries of the SW corner cell.
for i = 1:nih
    for j = 1:njh
        
        % this enforces MOM's NE convention
        [iu,ju] = deal(i+1,j);
        [iv,jv] = deal(i,j+1);
        
        %------------------------------------ matrix of full eddy flux
        F = [squeeze(fu(iu,ju,:))';
            squeeze(fv(iv,jv,:))'];
        F = reshape(F,[2*ntr,1]); % column vec [Fx1;Fy1;Fx2;..]
        
        %------------------------------------ matrix of tracer gradient
        % dCdx on u-grid, boundary points are set to ZERO!
        Gxu = [squeeze(dcdx_u(iu,ju,:))';
            zeros(1,ntr)           ];
        Gxu = reshape(Gxu,[2*ntr,1]);
        % dCdy on u-grid
        Gyu = [squeeze(dcdy_u(iu,ju,:))';
            zeros(1,ntr)           ];
        Gyu = reshape(Gyu,[2*ntr,1]);
        % dCdx on v-grid
        Gxv = [zeros(1,ntr);
            squeeze(dcdx_v(iv,jv,:))'];
        Gxv = reshape(Gxv,[2*ntr,1]);
        % dCdy on v-grid
        Gyv = [zeros(1,ntr);
            squeeze(dcdy_v(iv,jv,:))'];
        Gyv = reshape(Gyv,[2*ntr,1]);
        % form a 2*ndist-by-4 matrix
        G = [Gxu,Gyu,Gxv,Gyv];
        
        %------------------------------------ inverse
        % - G * K = F 
        if all(~isnan(G),'all') && all(~isnan(F),'all')
            
            % ~flg_ang or  rank(G,norm(G)/1e2) == 4
            % smaller the TOL, easier to full rank
%             if rank(G) == 4 % rank(G,norm(G)/1e2) == 4 or '/5e1'?
%                 K = - lsqminnorm(G,F);  % G\F different from pinv(G)*F !!
%             else
%                 K(1:4) = NaN;
%             end
            %tol = norm(G)/1e2;
    	    % K = - G \ F;
	    K = - pinv(G)*F;
            %K = - lsqminnorm(G,F);  % G\F different from pinv(G)*F !!
            Kxx(iu,ju) = K(1); Kxy(iu,ju) = K(2);
            Kyx(iv,jv) = K(3); Kyy(iv,jv) = K(4);
            
            %----------------- calc angles btw trac grad
%             angu(j,i) = angle_vectors(G(1,1:2), G(3,1:2));
%             angv(j,i) = angle_vectors(G(2,3:4), G(4,3:4));
            
        end
    end
end
