function [u,v] = p2uv_CG(u_p,v_p)
% 
% Interpolate p-grid u/v velocity fields onto u-/v-grids.
% 
[nih, njh] = size(u_p);
[niu, nju] = deal(nih+1, njh);
[niv, njv] = deal(nih, njh+1);

[u, v] = deal(zeros(niu,nju), zeros(niv,njv));

u(2:end-1,:) = ( u_p(1:end-1,:) + u_p(2:end,:) ) / 2;
v(:,2:end-1) = ( v_p(:,1:end-1) + v_p(:,2:end) ) / 2;

end