function [ x_bn ] = BatchNormalize( x,dim )

x_bn = DivideSTD(ReduceMean(x,dim),dim);

end
