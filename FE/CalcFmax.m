function [fmax,max_bin,X] = CalcFmax(x,fs,k,prop)

x = x - mean(x,1);
if strcmp(prop,'full')
    dims = 1:4;
else
    dims = 4;
end

fftsize = 512;
freq_vec = linspace(0,fs/2,fftsize/2);

X = fft(x(:,dims,:),fftsize);
X = X(1:fftsize/2,:,:);
[~,max_bin] = maxk(X,k,1);
fmax  = squeeze(freq_vec(max_bin));

if size(fmax,1) < size(fmax,2) 
    fmax = fmax';
end
end