function f = mismatch_func_2vec(x, fld1u, fld1v, fld2u, fld2v)
% fld1u/v, fld2u/v: x-/y- comp of the two vectors (same dim)
% x: a scalar, 
% f = integral of ||vec1 - x*vec2||^2
% 
if ndims(fld1u) ~= ndims(fld1v) || ndims(fld1u) ~= ndims(fld2u)
    error('Two input fields must have the same size!')
end
if ~isscalar(x)
    error('The ratio x should be a scalar for now!')
end

diffu = fld1u - x .* fld2u; 
diffv = fld1v - x .* fld2v; 
f_obj = (diffu + diffv).^2;
f = sum(f_obj, 'all', 'omitnan');

end
