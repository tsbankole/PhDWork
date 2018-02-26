function [t2,y,f_ov_ss, Rvss,u] = vandevusse_Nonlinear_separator(tspan,f_OV,y0)
Cas= 3;
Cbs = 1.117;
f_ov_ss =4/7;
interval = length(f_OV);

recov_ca = 0.9;
recov_cb = 0.1;
molar_vol = 0.1;

[t2,y]=ode45(@nonlinear_separator,tspan,y0,[],interval, tspan(2),f_OV);

for k = 1:interval
    u(find(t2>=(k-1)*tspan(2)/interval,1,'first'):find(t2<k*tspan(2)/interval,1,'last'))= f_OV(k);
end
u(end+1) = f_OV(end);

u = u'+f_ov_ss;
Rvss = f_ov_ss*molar_vol*(recov_ca*Cas + recov_cb*Cbs)/(1-molar_vol*(recov_ca*Cas + recov_cb*Cbs));

A = u.*molar_vol.*(recov_ca.*y(:,1) +recov_cb*y(:,2));
B = (1-molar_vol*(recov_ca*y(:,1) + recov_cb*y(:,2)));

Rv = A./B;
u = [u  Rv];

    function [xdot, Rv]=nonlinear_separator(t,x,interval, tend, f_ov)
        k1=5/6;
        k2=5/3;
        k3=1/6;
        f_ov_ss =4/7; %steady state overhead feed to volume ratio
        
        % extra parameters for separator
        vol_ratio = 1; % ratio of reactor to separator
        
        ca=x(1);
        cb=x(2);
        car = x(3);
        cbr = x(4);
        
        for i = 1: interval
            if t>=((i-1)/interval)*tend && t<i/interval *tend
                f_ov_u = f_ov(i);
            end
            if t>=tend
                f_ov_u = f_ov(end);
            end
        end
        
        fv = f_ov_u + f_ov_ss;
        Rv = fv*molar_vol*((recov_ca * ca)+(recov_cb * cb))/(1-molar_vol*((recov_ca * ca)+(recov_cb * cb)));
        
        caf=10;
        
        xdot(1)=fv*(caf-ca)-k1*ca-k3*ca^2 + Rv*(car-ca);
        xdot(2)=-fv*cb+k1*ca-k2*cb + Rv*(cbr-cb);
        xdot(3)= vol_ratio*(fv +Rv)*(ca -car) ;
        xdot(4)= vol_ratio*(fv +Rv)*(cb-cbr) ;
        %         xdot(5)= x(5) - fv*molar_vol*((recov_ca * ca)+(recov_cb * cb))/(1-molar_vol*((recov_ca * ca)+(recov_cb * cb)));
        %         xdot(6)= x(6) - (fv+Rv);
        
        xdot=xdot';
        
    end
end
