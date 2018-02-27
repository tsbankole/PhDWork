function [data1, u1, t1] = generate_data()
% to generate data for dcm
clear all
close all
Cas= 3;
Cbs = 1.117;
f_ov_ss =4/7;
Caf_ss=10;
tspan = [0 10];
interval = 2;
Ca0 = 2;
Cb0 = 1.117;
y0=[Ca0;Cb0;Cas;Cbs];
f_OV = [0 0.5];
[t2,y,f_ov_ss, Rvss,u] = vandevusse_Nonlinear_separator(tspan,f_OV,y0);
t_sp = [tspan(1):0.05:tspan(2)]';
y_sp = interp1(t2,y,t_sp,'spline') - repmat([Cas Cbs Cas Cbs], [length(t_sp) 1]);
u1 = interp1(t2,u,t_sp,'next') - repmat([f_ov_ss, Rvss], [length(t_sp) 1]);
u_sp = u_sp';
data1 = y_sp';
t1 = t_sp';
end
