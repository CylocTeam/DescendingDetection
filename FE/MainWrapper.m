utils = '';
DB_path = '';
close all
addpath(genpath(DB_path))
addpath(utils)

holder_vec      = {'A','D','O'};
building_id_vec = {'T','S','D'};
bundle_vec      = {'Bundle1','Bundle2','BundleX'};
device_vec      = {'D1','D2','D3'};
position_vec    = {'Pocket','Hand','Bag'};

fs_new = 15;
win_size = 3 * fs_new;
step_size = fix(0.2 * fs_new);

damping = 0.5;
fres = 3;
batch_norm = 1;
mefilt_labels = 0;
efficiency_bbox = [-5,95 ; 0,90 ; 5,95 ; 5,100 ; -5,100]';

features_struct = struct();
features_struct.    std.            is_used   = 1;    features_struct.std.args           = {   'norm'};
features_struct.    fmax.           is_used   = 1;    features_struct.fmax.args          = {1, 'norm'};
features_struct.    mad0.           is_used   = 1;    features_struct.mad0.args          = {   'norm'};
features_struct.    mad1.           is_used   = 1;    features_struct.mad1.args          = {   'norm'};
features_struct.    iqr.            is_used   = 1;    features_struct.iqr.args           = {   'norm'};
features_struct.    bandwidth.      is_used   = 1;    features_struct.bandwidth.args     = {   'norm'};
features_struct.    prctile5.       is_used   = 1;    features_struct.prctile5.args      = {5  'norm'};
features_struct.    prctile95.      is_used   = 1;    features_struct.prctile95.args     = {95,'norm'};
features_struct.    max_wavelet.    is_used   = 1;    features_struct.max_wavelet.args   = {1  'norm'};
features_struct.    sma.            is_used   = 1;    features_struct.sma.args           = {   'norm'};
features_struct.    thd.            is_used   = 1;    features_struct.thd.args           = {   'norm'};
features_struct.    pca.            is_used   = 1;    features_struct.pca.args           = {   'norm'};
features_struct.    trend.          is_used   = 1;    features_struct.trend.args         = {   'norm'};
features_struct.    mean.           is_used   = 0;    features_struct.mean.args          = {   'norm'};

%%
% clear result_struct
result_cell = {};
Smap = [];
Pmap = [];
for ibuilding=1:length(building_id_vec)
    building_id = building_id_vec{ibuilding};
    
    for iholder=1:length(holder_vec)
        holder = holder_vec{iholder};
        
        for ibundle=1:length(bundle_vec)
            bundle = bundle_vec{ibundle};
            
            for idevice=1:length(device_vec)
                device = device_vec{idevice};
                
                for iposition=1:length(position_vec)
                    position = position_vec{iposition};
                    
                    try
                        % close all
                        Main;  
                            result_cell(end+1,:) = { building_id , holder, bundle, device, position, recall, far, accuracy,error_imbalance, sill_pctile};
                            Smap = [ Smap ; p_egs];  
                            Pmap = [ Pmap ; p_map];
                    catch
                        disp({ building_id , holder, bundle, device, position})
                        continue
                    end
                end
            end
        end
    end
    disp(ibuilding)
end
% clr = colormap('jet');
% clr = clr(fix(linspace(1,size(clr,1),length(holder_vec))),:);
% clr_by_criteria = clr(findgroups(result_tbl.holder),:);


result_tbl = cell2table(result_cell,'VariableNames',{'building','holder','bundle','device','position','recall','false_alarm','acuuracy','error_imbalance','s25'});

accum_tbl = [accum_tbl , struct('features',features_struct,...
                                'result_table',result_tbl,...
                                'pctg590',sum(result_tbl.recall > 90 & result_tbl.false_alarm < 5) / size(result_tbl,1) * 100,...
                                'top_pctg', sum(inpolygon(result_tbl.error_imbalance,100*result_tbl.acuuracy,...
                                                          efficiency_bbox(1,:),efficiency_bbox(2,:)))/size(result_tbl,1)*100) ];
    
% calc ROC for s25
is_efficient = inpolygon(result_tbl.error_imbalance,100*result_tbl.acuuracy,...
                         efficiency_bbox(1,:),efficiency_bbox(2,:));
top_pctg = sum(is_efficient) / length(is_efficient);                    
[tpr,fpr,thr] = roc(is_efficient',result_tbl.s25');  
s25_recall_0fa = 100*tpr(find(fpr==0,1,'last'));
min_fa_thr = thr(find(fpr==0,1,'last'));

% calc threshold margin - 

% figure(1)
% scatter3(result_tbl.false_alarm, result_tbl.recall,result_tbl.s25,[],result_tbl.s25,'filled')
% colormap('jet')
% grid minor
% xlabel('False Alarm [%]')
% ylabel('Recall [%]')
% feature_list = features_tbl.Properties.VariableNames;
% parsed_labels = cellfun(@(x) [x,', '],feature_list,'UniformOutput',false);
% parsed_labels = strrep(parsed_labels,'_','__');
% 
% title(sprintf(['Features: ' [parsed_labels{:}],'\n Efficiency: %.2f [%%]'],accum_tbl(end).pctg590))

figure(2)
scatter3(result_tbl.error_imbalance, 100*result_tbl.acuuracy,result_tbl.s25,[],result_tbl.s25,'filled')
hold all
colormap('jet')
grid minor
xlabel('Error Imbalance [%]')
ylabel('Accuracy [%]')
feature_list = features_tbl.Properties.VariableNames;
parsed_labels = cellfun(@(x) [x,', '],feature_list,'UniformOutput',false);
parsed_labels = strrep(parsed_labels,'_','__');

title(sprintf(['Features: ' [parsed_labels{:}],'\n Efficiency: %.2f [%%]\n s25 Recall: %.2f [%%] @ TH=%.2f\n Total Efficiency: %.2f [%%]'],...
              100*top_pctg,s25_recall_0fa,min_fa_thr,top_pctg*s25_recall_0fa))
plot([efficiency_bbox(1,:),efficiency_bbox(1,1)],[efficiency_bbox(2,:),efficiency_bbox(2,1)],'r','LineWidth',3)
colorbar    



%% 

figure; imagesc('XData',s_th_vec,'CData',Smap);%(Smap(:,end)<0.5,:))

clear profitabillity_pctg_unbiased
for icol=1:size(Smap,2)
    curr_s_unb = Smap(:,icol);
    profitabillity_pctg_unbiased(icol) = mean(curr_s_unb(~isnan(curr_s_unb)) > 0.5 ); % unbiased
    profitabillity(icol)               = mean(curr_s_unb(~isnan(curr_s_unb))) - 0.5; 
end

figure; plot(s_th_vec,100*profitabillity_pctg_unbiased,'LineWidth',2)
hold all
        yyaxis right
        plot(s_th_vec, 100*profitabillity,'LineWidth',2)
grid minor
grid on


%% 
figure; imagesc('XData',s_th_vec,'CData',Pmap);%(Smap(:,end)<0.5,:))

clear profitabillity_pctg_unbiased
for icol=1:size(Pmap,2)
    curr_s_unb = Pmap(:,icol);
    profitabillity_pctg_unbiased(icol) = mean(curr_s_unb(~isnan(curr_s_unb)) > 0.5 ); % unbiased
    profitabillity(icol)               = mean(curr_s_unb(~isnan(curr_s_unb))) - 0.5; 
end

figure; plot(s_th_vec,100*profitabillity_pctg_unbiased,'LineWidth',2)
hold all
        yyaxis right
        plot(s_th_vec, 100*profitabillity,'LineWidth',2)
grid minor
grid on
