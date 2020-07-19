function [ dims ] = ParseInputDimentions( spec )

if isnumeric(spec)
    dims = spec;
    return
end

switch spec
    case 'full'
        dims = 1:6;
    case 'norm'
        dims = 4;
    case 'cart'
        dims = 1:3;
    case 'cart-norm'
        dims = 1:4;
    case 'azimuth'
        dims = 5;
    case 'elevation'
        dims = 6;
    case 'sphere'
        dims = 4:6;
    otherwise
        dims = spec;
end

end

