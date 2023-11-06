function [u_p,v_p] = uv2p_CG(u,v)
% 
%  u: [niu-nju]
%  v: [niv,njv]
% 
% [niu, nju] = size(u);
% [nip, njp] = deal(niu-1, nju);
% [u_p,v_p] = deal(NaN * zeros(nip, njp);
u(isnan(u)) = 0;
v(isnan(v)) = 0;

u_p = (u(1:end-1,:) + u(2:end,:)) / 2;
v_p = (v(:,1:end-1) + v(:,2:end)) / 2;

end
