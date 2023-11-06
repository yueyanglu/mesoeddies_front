function f = mismatch_func(x, fld1, fld2)
% fld1, fld2: matrices of the same dim
% x: a scalar, 
% f = integral of (fld1 - x*fld2) 
% 
if ndims(fld1) ~= ndims(fld2)
    error('Two input fields must have the same size!')
end
if ~isscalar(x)
    error('The ratio x should be a scalar for now!')
end

f_obj = (fld1 - x .* fld2).^2;
f = sum(f_obj, 'all', 'omitnan');

end
