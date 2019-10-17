function result=is_finite_number(foo_array)
result=prod(isfinite(foo_array) & isnumeric(foo_array))==1;