function [data_enframed,labels_enframed] = EnframeData(acc_data,labels,win_size,step_size)
[N,D] = size(acc_data);
NumFrames = fix(N / step_size - 1);

clear data_enframed
for d=1:D
    data_enframed(:,d,:) = buffer(acc_data(:,d),win_size,win_size-step_size,'nodelay');
end

labels_enframed = buffer(labels,win_size,win_size-step_size,'nodelay');
labels_enframed = mode(labels_enframed,1);
labels_enframed = labels_enframed(:);
end

% figure; plot(labels_enframed)
% figure; plot(labels)
