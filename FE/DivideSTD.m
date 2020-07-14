function [ x_divided ] = DivideSTD(x,dim)

expand_size = repmat({1},size(x)); 
expand_size{dim} = size(x,dim);

x_divided = x ./ repmat(std(x,[],dim),expand_size{:});

end
