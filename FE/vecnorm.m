function [ nvorm ] = vecnorm( v,p,dim )
nvorm = (sum(v.^p,dim)).^(1/p);
end

