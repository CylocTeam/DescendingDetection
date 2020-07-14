function [] = PlotVideoScatter(acc_data,features,labels)
clrs = {'b','r','g'};
Npts = size(features,1);
acc_norm = squeeze(vecnorm(acc_data,2,2));
fs = 15;
close all
figure(1)

grid minor
figure(2)
grid minor

% last_down = find(labels==1,1,'last')+ 1;

% for ipt=last_down:Npts
for ipt=1:Npts

    curr_acc = acc_norm(:,ipt);
    curr_pt = features(ipt,:);
    
    figure(1);
    cla
    plot(curr_acc,clrs{labels(ipt)+1})
    ylim([-5 35])
    xlim([-20 65])
    
    figure(2)
    scatter(curr_pt.std_norm,curr_pt.fmax_norm,'o','filled','MarkerFaceColor',clrs{labels(ipt)+1}); hold all
    hold all
    xlim([-3.5 3.5])
    ylim([-3.5 3.5])
    
    
    figure(3)
    plot(linspace(-15/2,15/2,512),pow2db(abs(fftshift(fft(curr_acc-mean(curr_acc),512)))),clrs{labels(ipt)+1})
    ylim([-20 30])
    
%     figure(4)
%     fft_vec = pow2db(abs(fft(curr_acc-mean(curr_acc),512)));
%     fft_vec = (fft_vec - mean(fft_vec))/std(fft_vec);
%     [f,x] = ecdf(fft_vec(2:end/2));
%     plot(x,f,clrs{labels(ipt)+1}); hold on;
    %     pause(0.2)
%     thd_val = thd(curr_acc,fs);
%     figure(4)
%     plot(ipt,thd_val,'o','Color',clrs{labels(ipt)+1}); 
%     hold on;
    
%     figure(4);  
%     plot(ipt,mad(cumsum(curr_acc)-mean(cumsum(curr_acc)),0),'o','Color',clrs{labels(ipt)+1})
%     hold on
%     axis tight
end

figure(2)
xlabel('STD')
ylabel('FMAX')
end

