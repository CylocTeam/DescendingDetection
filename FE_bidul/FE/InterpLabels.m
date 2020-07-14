function [labels_interp,t_interp] = InterpLabels(labels, t, fs_new)

t_strat = t(1);
t_end = t(end);
t_interp = t_strat : 1/fs_new : t_end;
labels_interp = round(interp1(t,labels,t_interp,'spline'));
labels_interp(labels_interp < 0) = 0;
labels_interp(labels_interp > 2) = 2;

end