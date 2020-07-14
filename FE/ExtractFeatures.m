function [features] = ExtractFeatures(acc_data,fs,feature_usage_struct,batch_norm)

feature_bank.std         = @(x,prop)     CalcSTD(x,prop);
feature_bank.fmax        = @(x,nf,prop)  CalcFmax(x,fs,nf,prop);
feature_bank.prctile5    = @(x,pct,prop) CalcPctile(x,5,prop);
feature_bank.prctile95   = @(x,pct,prop) CalcPctile(x,95,prop);
feature_bank.iqr         = @(x,prop)     CalcIQR(x,prop);
feature_bank.mad0        = @(x,prop)     CalcMad(x,0,prop);
feature_bank.mad1        = @(x,prop)     CalcMad(x,1,prop);
feature_bank.max_wavelet = @(x,nw,prop)  CalcWaveletCoeffs(x,nw,prop);
feature_bank.bandwidth   = @(x,prop)     CalcBandWidth(x,fs,prop);
feature_bank.sma         = @(x,prop)     CalcSMA(x,prop);
feature_bank.thd         = @(x,prop)     CalcTHD(x,fs,prop);
feature_bank.pca         = @(x,prop)     CalcPCA(x,prop);
feature_bank.mean        = @(x,prop)     CalcMean(x,prop);
feature_bank.trend       = @(x,prop)     CalcTrend(x,fs,prop);

func_bank = struct2cell(feature_bank);
func_names = fieldnames(feature_bank);

%% normalize data
% calc |a|
% acc_norm = vecnorm(acc_data,2,2);
acc_norm = vecnorm(acc_data,2,2);

acc_data_ext = cat ( 2, acc_data , acc_norm );

%% calc features
fnames = fieldnames( feature_usage_struct );
lgd_types = {'x','y','z','norm'};
features = [];
all_headers = {};

for ifunc=1:length(fnames)
    curr_name = fnames{ifunc};
    curr_feature_struct = feature_usage_struct.(curr_name);
    
    % check usage flag & validity
    if ~curr_feature_struct.is_used
        continue
    end
    
    if strcmp(curr_feature_struct.args{end},'full')
        feature_headers = cellfun(@(t) [curr_name,'_',t],lgd_types,'UniformOutput',false);
    else
        feature_headers = [curr_name,'_',lgd_types{4}];
    end
    
    % calc_features
    curr_features =  feature_bank.(curr_name)(acc_data_ext,curr_feature_struct.args{:});
    
    % concat headers
    features     = horzcat(features, curr_features);
    all_headers  = horzcat(all_headers,feature_headers);
end

features(end,:) = features(end-1,:) ;   % because of 0 padding at buffer()


%%
if batch_norm
%     features = (features - mean(features,1))./ std(features,[],1);
    features = BatchNormalize(features,1);
end

% std_idx = contains(all_headers,'std_norm');
std_idx = ismember(all_headers,'std_norm');

if any(std_idx)  % actualy should by exacly single "1"
    min_std = prctile(features(:,std_idx),20);         % assume lower 10% of std is stay
    features(features(:,std_idx) < min_std ,:) = 0;
%     features.std(features.std < min_std ) = 0;
end

if batch_norm  % again
%     features = (features - mean(features,1))./ std(features,[],1);
      features = BatchNormalize(features,1);
end

features =  array2table(features,'VariableNames',all_headers);

end

function [std_val]   = CalcSTD(x,prop)

if strcmp(prop,'full')
    dims = 1:4;
else    % norm
    dims = 4;
end
std_val = squeeze(std(x(:,dims,:),[],1));
if size(std_val,1) < size(std_val,2)
    std_val = std_val';
end
end
function [mad_val]   = CalcMad(x,flag,prop)

if strcmp(prop,'full')
    dims = 1:4;
else    % norm
    dims = 4;
end
mad_val = squeeze(mad(x(:,dims,:),flag,1));
if size(mad_val,1) < size(mad_val,2)
    mad_val = mad_val';
end
end
function [pctl]      = CalcPctile(x,pct,prop)
% x = x - mean(x,1);
x = ReduceMean(x,1);
if strcmp(prop,'full')
    dims = 1:4;
else    % norm
    dims = 4;
end
pctl = squeeze(prctile(x(:,dims,:),pct,1))';
if size(pctl,1) < size(pctl,2)
    pctl = pctl';
end
end
function [iqr]       = CalcIQR(x,prop)
if strcmp(prop,'full')
    dims = 1:4;
else    % norm
    dims = 4;
end
iqr = squeeze(prctile(x(:,dims,:),75,1) - prctile(x(:,dims,:),25,1));
if size(iqr,1) < size(iqr,2)
    iqr = iqr';
end
end
function [max_coeff] = CalcWaveletCoeffs(x,nw,prop)
% x = x - mean(x,1);
x = ReduceMean(x,1);

if strcmp(prop,'full')
    dims = 1:4;
else    % norm
    dims = 4;
end
[Lo_D,Hi_D] = wfilters('haar','d');

[N,M,L] = size(x);
Ncoeffs = fix(log2(N));
x = x(:,dims,:);

clear c
for icol=1:size(x,2)
    for iwin=1:size(x,3)
        c(iwin,icol,:) = wavedec(x(:,icol,iwin), Ncoeffs,Lo_D,Hi_D);
    end
end
% [~,max_coeff] = maxk(c,nw,3);
[~,max_coeff] = max(c,[],3);

if size(max_coeff,1) < size(max_coeff,2)
    max_coeff = max_coeff';
end

end
function [bw]        = CalcBandWidth(x,fs,prop)
if strcmp(prop,'full')
    dims = 1:4;
else    % norm
    dims = 4;
end

[~,fmax_idx,X] = CalcFmax(x,fs,1,prop);
bw_drop = 3;
fftsize_single_sided = size(X,1);
freq_vec = linspace(0,fs/2,fftsize_single_sided);
bw = zeros(size(x,3),1);

for icol=1:size(X,2)
    
    for iwin=1:size(X,3)
        curr_fft = pow2db(abs(X(:,icol,iwin)));
        max_idx = fmax_idx(:,icol,iwin);
        
        % right drop
        right_idx = max_idx;
        drop = 0;
        while(drop < bw_drop && right_idx < fftsize_single_sided)
            right_idx = right_idx + 1;
            drop = curr_fft(max_idx) - curr_fft(right_idx) ;
        end
        
        % left drop
        left_idx = max_idx;
        drop = 0;
        while(drop < bw_drop && left_idx > 0)
            left_idx = left_idx - 1;
            drop = curr_fft(max_idx) - curr_fft(left_idx) ;
        end
        
        bw(iwin,icol) = freq_vec(right_idx) - freq_vec(left_idx);
    end
end

end
function [sma]       = CalcSMA(xyz,prop)

x = xyz(:,1:3,:);

N = size(x,1);
sma = squeeze(trapz(sum(abs(x),2),1))/N;
if size(sma,1) < size(sma,2)
    sma = sma';
end
sma(end) = sma(end-1);  % because of 0 padding of buffer

end
function [thd_x]     = CalcTHD(x,fs,prop)
% x = x - mean(x,1);
x = ReduceMean(x,1);

if strcmp(prop,'full')
    dims = 1:4;
else    % norm
    dims = 4;
end
x = x(:,dims,:);
clear thd_x
for icol=1:size(x,2)
    for iwin=1:size(x,3)
        thd_x(iwin,icol) = thd(x(:,icol,iwin),fs);
    end
end
thd_x(isinf(thd_x) | isnan(thd_x)) = -30;
thd_x = movvar(thd_x,9,1);
end
function [tsquared]  = CalcPCA(x,prop)
% x = x - mean(x,1);
x = ReduceMean(x,1);

if strcmp(prop,'full')
    dims = 1:4;
else    % norm
    dims = 4;
end
x = x(:,dims,:);
clear tsquared
for icol=1:size(x,2)
    [~,~,~,tsquared(:,icol)] = pca(squeeze(x(:,icol,:))');
    
end
end
function [mn]        = CalcMean(x,prop)
if strcmp(prop,'full')
    dims = 1:4;
else    % norm
    dims = 4;
end
mn = squeeze(mean(x(:,dims,:),1));
if size(mn,1) < size(mn,2)
    mn = mn';
end
end
function [trend_std]      = CalcTrend(x,fs,prop)
if strcmp(prop,'full')
    dims = 1:4;
else    % norm
    dims = 4;
end
x = x(:,dims,:);

clear mean_v
for iwin=1:size(x,3)
    for icol=1:size(x,2)
        curr_x = x(:,icol,iwin);
        detrend_x = detrend(curr_x);
        
        curr_v = cumtrapz(0:1/fs:(length(detrend_x)-1)/fs,detrend_x);
        mean_v(iwin,icol) = mean(curr_v);
        
    end
end
trend_std = movvar(mean_v,9,1);
end


