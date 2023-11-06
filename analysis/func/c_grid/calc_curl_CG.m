function curlF = calc_curl_CG(Fx,Fy,dxBu,dyBu)
%Calculate the horizontal curl of a vector F = (Fx, Fy), which is defined
%   as (z X del) \cdot F = (-d/dy, d/dx) \cdot F = - dFu/dy + dFv/dx 
% 
% 
% INPUT:
% Fx - [lonq,lath]
% Fy - [lonh,latq]
% dxBu - [lonq,latq]
% dyBu - [lonq,latq]
% 
% OUTPUT: 
% (horizontal) curl of F [lonq,latq]
% 

[niq, njq] = size(dxBu);
% [njq2,niq2] = size(scqx);
% if njq2 ~= njq || niq2 ~= niq
%     error('Grid and flux dimensions do not match!');
% end

curlF = NaN * zeros(niq,njq);

[j_op, i_op] = deal(2:njq-1, 2:niq-1);

% [niq-2, njq-2]
dvdx = (Fy(2:end,2:end-1) - Fy(1:end-1,2:end-1)) ./ dxBu(i_op,j_op);
dudy = (Fx(2:end-1,2:end) - Fx(2:end-1,1:end-1)) ./ dyBu(i_op,j_op);


curlF(i_op,j_op) = dvdx - dudy;