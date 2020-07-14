function [ x_reduced ] = ReduceMean(x,dim)

expand_size = repmat({1},size(x)); 
expand_size{dim} = size(x,dim);

x_reduced = x - repmat(mean(x,dim),expand_size{:});

end

