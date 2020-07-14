function [acc_interp,t_interp] = InterpSignal(acc, t, fs_new)

t_strat = t(1);
t_end = t(end);
t_interp = t_strat : 1/fs_new : t_end;
acc_interp = interp1(t,acc,t_interp,'spline');

end