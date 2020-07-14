%% 0 - Read Data + validate uniquness
featuers_file = dir(fullfile(DB_path,holder, building_id, bundle, device,sprintf('%s_%s_%s_%s_%s.csv',holder,building_id,bundle,device,position)));
raw_data = readtable(fullfile(featuers_file.folder,featuers_file.name));
raw_data.time = raw_data.time - raw_data.time(1);  % [sec] 

[~,unq_time] = unique(raw_data.time);
raw_data = raw_data(unq_time,:);

%% 1 LPF + Resample
lpf_data = TwoPoleLPF(raw_data.time, [raw_data.accx,raw_data.accy, raw_data.accz],fres,damping);
[lpf_resampled_data,t_interp] = InterpSignal(lpf_data,raw_data.time,fs_new);
[lpf_resampled_label] = InterpLabels(raw_data.label,raw_data.time,fs_new);

%% 2 - Enframe
[data_enframed,labels_enframed] = EnframeData(lpf_resampled_data,lpf_resampled_label,win_size,step_size);

%% 3 - Feature Extraction
features_tbl = ExtractFeatures(data_enframed,fs_new,features_struct,batch_norm);
% features_tbl.gsp =  vecnorm(gradient(table2array(features_tbl)),2,2).*features_tbl.std_norm;
% features_tbl.gsp = (features_tbl.gsp - mean(features_tbl.gsp)) / std(features_tbl.gsp);
%% 4 - Cluster
lbls_kmeans = (double(labels_enframed > 1)) ;  % [binary label]


ms_features = MeanShift(table2array(features_tbl),1,1,1,1);
class_kmeans = kmeans(ms_features,2,'Replicates',10,'Distance','cityblock');

% class_kmeans = kmeans(table2array(features_tbl),2,'Replicates',10,'Distance','cityblock');
% gmmodel = fitgmdist([table2array(features_tbl),lbls_kmeans],2);

%%
% [~,max_avg_norm_idx] = max(splitapply(@(x) mean(vecnorm(x,1,2)),table2array(features_tbl),class_kmeans));
% is_down_class = (class_kmeans == max_avg_norm_idx);  % assume 1st point has label "0"

% [Fx,Fy] = gradient(table2array(features_tbl));
% figure; 
% % plot(10*gradient(features_tbl.std_norm).*features_tbl.std_norm)
% % plot(10*movvar(gradient(features_tbl.std_norm).*features_tbl.std_norm,25))
% hold all
% plot(labels_enframed)
% 
% figure; 
% plot(vecnorm(Fx,2,2).*features_tbl.std_norm)
% hold all
% plot(features_tbl.std_norm)
% plot(labels_enframed)
% 
% 
% figure; 
% plot(vecnorm(Fx,2,2))
% hold all
% plot(labels_enframed)


% post proccess - rename classes
% [~,max_std_idx] = max(features_tbl.std_norm);
% is_down_class = (class_kmeans == class_kmeans(max_std_idx)); 
% if is_down_class(1)   
%     is_down_class = ~is_down_class;
% end
is_down_class = (class_kmeans ~= class_kmeans(1));  % assume 1st point has label "0"
named_class = double(is_down_class);

%% filter clusters over time
if mefilt_labels
    named_class = medfilt1(named_class,...
                           floor( 4.5 / (step_size/fs_new) / 2)*2 + 1);  % closest odd number of x is:  1 + 2*floor(x/2)
end
[s] = silhouette(table2array(features_tbl),class_kmeans,'cityblock');
sill_pctile = prctile(s,25);
% ecdf(s)

%% Calc Accuracy
pdist_lbl = (lbls_kmeans == lbls_kmeans');
pdist_kmeans = (named_class == named_class');
N = length(lbls_kmeans);
paired_accuracy = sum(sum(triu( pdist_lbl == pdist_kmeans))) / ((N+1)*N/2) ; 
% accuracy = sum(named_class == lbls_kmeans) / length(lbls_kmeans);
% fprintf('\n\nAcc: %.2f\nPaired-Acc: %.2f\n\n',accuracy, paired_accuracy)

[miss_accuracy,cm] = confusion(lbls_kmeans',named_class' );
% disp(cm)
accuracy = 1 - miss_accuracy;
recall = cm(2,2) / sum(cm(2,:)) * 100 ;
far    = cm(1,2) / sum(cm(1,:)) * 100 ;
error_imbalance = ( cm(2,1) - cm(1,2) ) / length(lbls_kmeans) * 100;

%% Plot Predictions
% figure;
% % plot(t_interp,squeeze(vecnorm(data_enframed,2,2))); hold all;
% plot(lbls_kmeans,'.','MarkerSize',20); hold all
% plot(named_class+0.25,'.','MarkerSize',20); 
% plot(find(named_class ~= lbls_kmeans),named_class(named_class ~= lbls_kmeans)+0.25,'ok','MarkerSize',10)
% plot(s,'.k')
% ylim([-2 5])
% grid minor
% i=1;

% figure; 
% % ksdensity(s(named_class == lbls_kmeans),'Bandwidth',0.02,'Kernel','epanechnikov')  % accuracy
% % hold all
% % ksdensity(s(named_class==0 & lbls_kmeans==1),'Bandwidth',0.02,'Kernel','epanechnikov')  % MD
% % ksdensity(s(named_class==1 & lbls_kmeans==0),'Bandwidth',0.02,'Kernel','epanechnikov')  % FA
% 
% ecdf(s(named_class==1 == lbls_kmeans==1)); % TP
% hold all
% ecdf(s(named_class==0 & lbls_kmeans==1));  % FN
% ecdf(s(named_class==0 & lbls_kmeans==0));  % TN
% ecdf(s(named_class==1 & lbls_kmeans==0));  % FP
% legend({'TP','MD','TN','FA'})

% calc bayesian error prob. using s
s_th_vec = -0.5:0.005:1;
for iis = 1:length(s_th_vec)
    s_th = s_th_vec(iis);
    p_s = sum(s < s_th)/length(s);
    p_e = miss_accuracy;
    p_sge = sum(s(named_class ~= lbls_kmeans) < s_th ) / sum(named_class ~= lbls_kmeans);
    p_egs(iis) = p_sge * p_e / p_s;
    p_map(iis) = p_sge * p_e ;  
end
% plot(s_th_vec, p_egs); hold all;


%% Plot Scatter for 2 features
% figure;
% for ilabel=0:2
%     label_idxs = (labels_enframed == ilabel);
%     try
%         scatter3(features_tbl.std_norm(label_idxs),features_tbl.fmax_norm(label_idxs),features_tbl.pca_norm(label_idxs),'o','filled'); hold all
%     catch
%         scatter (features_tbl.std_norm(label_idxs),features_tbl.fmax_norm(label_idxs),'o','filled'); hold all
%     end
% end
% xlabel('STD')
% ylabel('FMAX')
% zlabel('MAX')
% grid minor
% 
% figure;
% for ilabel=1:3
%     label_idxs = (class_kmeans == ilabel);
%     try
%         scatter3(features_tbl.std_norm(label_idxs),features_tbl.fmax_norm(label_idxs),features_tbl.pca_norm(label_idxs),'o','filled'); hold all
%     catch
%         scatter (features_tbl.std_norm(label_idxs),features_tbl.fmax_norm(label_idxs),'o','filled'); hold all
%     end
% end
% xlabel('STD')
% ylabel('FMAX')
% zlabel('MAX')
% grid minor
% i=1;
%% video scatter
% PlotVideoScatter(data_enframed,features_tbl,labels_enframed)

%% Features Vs Time
% PlotFeaturesVsTime(features_tbl,labels_enframed)
% figure; plot(labels_enframed); hold all
% plot(features_tbl.trend_norm); 
% plot(conv(hamming(7)/sum(hamming(7)),features_tbl.trend_norm),'LineWidth',2)
% plot(TwoPoleLPF(1:step_size/fs_new:size(features_tbl,1)*step_size/fs_new,features_tbl.trend_norm,0.5,damping),'LineWidth',2)
% 
% grid minor



