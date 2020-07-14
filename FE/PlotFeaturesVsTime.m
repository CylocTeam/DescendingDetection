function [] = PlotFeaturesVsTime(features,labels)

figure;
plot(table2array(features));
hold all
plot(10*(labels-1),'k')
grid on
xlabel('sample')
ylabel('features/label')
legend(horzcat(features.Properties.VariableNames,'label'))


% figure;
% grad_vec = diff([0, 0 ; features.std, features.fmax1]);
% [grad_theta,grad_rho] = cart2pol(grad_vec(:,1),grad_vec(:,2));
%  
% plot([grad_theta,grad_rho])
% hold all
% plot(10*(labels-1),'k')
% grid on
% xlabel('sample')
% ylabel('d(features)/dt , label')


end

