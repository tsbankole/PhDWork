clear
Cas= 3;
Cbs = 1.117;
f_ov_ss =4/7;
Caf_ss=10;
tspan = [0 10];
interval = 2;
Ca0 = 2;
Cb0 = 1.117;
y0=[Ca0;Cb0;Ca0;Cb0];
f_OV = [0 0.5];
[t2,y,f_ov_ss, Rvss,u] = vandevusse_Nonlinear_separator(tspan,f_OV,y0);
t_sp = [tspan(1):0.05:tspan(2)]';
y_sp = interp1(t2,y,t_sp,'spline');
y_sp = y_sp' - repmat([Cas; Cbs; Cas; Cbs], [1 size(y_sp,1)]);
t_sp = t_sp';
clc
