function [] = JPCtestcase1_Bankole()
% in the paper Exploiting Connectivity Structures for Decomposing Process Plants, the van de vusse reactor was approximated with a bilinear/second order nonlinear model and 
% the ensuing parameters A,B,H which form the parameter set are estimated and used for identifying structural connectivity. This MATLAB code reproduces the figures for which noisy data (avg of 30 here)
% is used for estimation. Plots are also provided to show how the parameter vector differs from the true vector and also the dataset variation with the underlying signal

pack = generate_data();
param(6).c  = pack.d;
param(10).c  = param(6).c;
fac = [0.0200    0.0050    0.0200    0.0200];
for ii = 1:4
    param(6).c(ii,2:end) = pack.d(ii,2:end) + fac(ii)*randn(size(pack.d(ii,2:end)));
end

pack.d=param(6).c ;

%%
toldiff = 0.005;
iter=1;
param(5).c =0.5;
param(8).c = false;
packD = dcm(pack);
Guess(:,1) = [reshape(packD.p1(2:end,2:end)',[16 1]);reshape(packD.p2(2:end,:,1)',[20 1]);reshape(packD.p2(2:end,:,2)',[20 1]);reshape(packD.p3(2:end,2:end,1)',[16 1])];
curentdiff(iter) = norm(packD.d(2:end,:) - param(6).c )
pack.d=param(5).c *packD.d(2:end,:) + (1-param(5).c )*param(6).c ;
XQ = rand(5,length(pack.t),10);
XQ(:,:,1)=packD.d;

while true
    iter=iter+1
    do = packD.d;
    packD = dcm(pack);
    Guess(:,iter) = [reshape(packD.p1(2:end,2:end)',[16 1]);reshape(packD.p2(2:end,:,1)',[20 1]);reshape(packD.p2(2:end,:,2)',[20 1]);reshape(packD.p3(2:end,2:end,1)',[16 1])];
    XQ(:,:,iter)=packD.d;
    curentdiff(iter) = norm(packD.d(2:end,:) - param(6).c );
    
    if curentdiff(iter)>curentdiff(iter-1)
        param(8).c = true;
        iter=iter-1;
        packD.d=do;
        param(5).c  = param(5).c -0.1;
    end
    
    if param(5).c <0.1 || param(5).c >0.9|| iter>25
        break
    end
    
    if ~param(8).c
        if curentdiff(iter-1)-curentdiff(iter)<toldiff
            param(5).c  = param(5).c +0.1;
        end
    end
    
    
    
    pack.d = param(5).c *packD.d(2:end,:) + (1-param(5).c )*param(6).c;
    
    param(8).c = false;
end


% 
% poj = dataraw-data1;
% for ij=1:4
% snri(ij) = snr(data1(ij,:),poj(ij,:));
% end
% disp('average signal to noise ratio = ')
% mean(snri) % establish signal to noise ratio

truevecparams = [-2.63000000000000;0;0.220000000000000;0;0.830000000000000;-2.46000000000000;0;0.220000000000000;0.790000000000000;0;-0.790000000000000;0;0;0.790000000000000;0;-0.790000000000000;7;-1;0;0;0;-1.11700000000000;0;-1;0;0;0;1;0;-1;0;0;0;1;0;-1;0;-1;0;1;0;0;0;-1;0;1;0;1;0;-1;0;0;0;1;0;-1;-0.166700000000000;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
for kk=1:size(Guess,2)
   normvec(kk)= norm(truevecparams-Guess(:,kk));
end
%
for ki= 1:2
figure(ki)
plot(pack.t, packD.d(ki+1,:),'r-', pack.t, param(6).c(ki,:),'*', pack.t,param(10).c(ki,:))
legend('estimate', 'noisy signal','true signal')
end

figure(3)
plot(normvec, 'b-*')
ylabel('$\|\theta_{true} - \hat{\theta}\|$','Interpreter','latex')

figure(4)
plot(curentdiff, 'b-*')
ylabel('$\|y - \hat{y}\|$','Interpreter','latex')


%
% set(0,'defaultlinelinewidth',1)
% set(0,'defaultaxeslinewidth',1)
% text = ['C_P   '; 'C_Q   ';'C_{Pr}';'C_{Qr}']
% set(gcf, 'Units','centimeters', 'Position',[5 5 9.5 9.2])
% for i=1:4
% figure(i)
% set(gcf, 'Units','centimeters', 'Position',[5 5 9.5 9.2])
% plot(t,dataraw(i,:),'k*',t,xq(i+1,:),'r-',t,data1(i,:),'b-.')
% legend('noisy signal','estimate','true signal','location','southeast')
% xlabel('t')
% ylabel(text(i,:))
% set(gcf,'color','none')
% set(gca, 'Fontsize',11)
% exportname = sprintf('noisyem%dsnr2.png', i);
% fname = sprintf('noisyem%dsnr2.fig', i);
% export_fig (exportname,'-painters','-r600')
% savefig(fname)
% end



    function pack  = generate_data()
        % to generate data for dcm
        clear all
        close all
        vdv.t = [0 10];
        vdv.ex = [0 0.5];
        vdv.in = [2; 1.117; 0; 0];
        vdvO = vandevusse_Nonlinear_separator(vdv);
        pack.t = [vdv.t(1):0.05:vdv.t(2)]';
        pack.u = transpose(interp1(vdvO.t, vdvO.u, pack.t, 'linear')   - repmat([vdvO.fss, vdvO.rss], [length(pack.t) 1]));
        pack.d = transpose(interp1(vdvO.t, vdvO.o, pack.t, 'spline') - repmat([3 vdv.in(2) 3 vdv.in(2)], [length(pack.t) 1]));
        pack.t = transpose(pack.t);
    end

    function vdvO = vandevusse_Nonlinear_separator(vdv)
        
        vdvO.fss =4/7;
        
        Inn.int = length(vdv.ex);
        Inn.lst = vdv.t(2);
        Inn.ex = vdv.ex;
        
        rca = 0.9;
        rcb = 0.1;
        mvl = 0.1;
        
        [vdvO.t,vdvO.o]=ode45(@nonlinear_separator,vdv.t,vdv.in,[],Inn);
        
        for k = 1:Inn.int
            vdvO.u(find(vdvO.t>=(k-1)*vdv.t(2)/Inn.int,1,'first'):find(vdvO.t<k*vdv.t(2)/Inn.int,1,'last'))= vdv.ex(k);
        end
        
        vdvO.u(end+1) = vdv.ex(end);
        vdvO.u = vdvO.u'+vdvO.fss;
        
        vdvO.rss = vdvO.fss*mvl*(rca*vdv.in(3) + rcb*vdv.in(4))/(1-mvl*(rca*vdv.in(3) + rcb*vdv.in(4)));
        
        vdvO.u = [vdvO.u  vdvO.u.*mvl.*(rca.*vdvO.o(:,1) +rcb*vdvO.o(:,2))./(1-mvl*(rca*vdvO.o(:,1) + rcb*vdvO.o(:,2)))];
        
        function [xdot, Rv]=nonlinear_separator(t,x,Inn)
            k1=5/6;
            k2=5/3;
            k3=1/6;
            
            vol_ratio = 1;
            
            ca=x(1);
            cb=x(2);
            car = x(3);
            cbr = x(4);
            
            for i = 1: Inn.int
                if t>=((i-1)/Inn.int)*Inn.lst && t<i/Inn.int*Inn.lst
                    f_ov_u = Inn.ex(i);
                end
                if t>=Inn.lst
                    f_ov_u = Inn.ex(end);
                end
            end
            
            fv = f_ov_u + vdvO.fss;
            Rv = fv*mvl*((rca * ca)+(rcb * cb))/(1-mvl*((rca * ca)+(rcb * cb)));
            
            caf=10;
            
            xdot(1)=fv*(caf-ca)-k1*ca-k3*ca^2 + Rv*(car-ca);
            xdot(2)=-fv*cb+k1*ca-k2*cb + Rv*(cbr-cb);
            xdot(3)= vol_ratio*(fv +Rv)*(ca -car) ;
            xdot(4)= vol_ratio*(fv +Rv)*(cb-cbr) ;
            
            xdot=xdot';
            
        end
    end



    function packD =  dcm(pack)
        %% Utilites
        runtime = zeros(20,1);
        d = pack.d;
        u = pack.u;
        %% ismember check
        % idx = arrayfun(@(x)find(A==x,1),B);
        setAA = 57;
        setB = sort([17 22 18 24 30 36 28 34 38 44 50 56 40 46 48 54]);
        %% define the set of nonzero parameters in set A/B
        delt = pack.t(2)-pack.t(1); % evenly spaced d points
        lengthu = size(u,1);
        N = size(d,1);
        N1 = N+1;
        N2 = N^2;
        
        if ~exist('setA','var')
            setA = 1:N2;
            %     flag = true;
        end
        
        numsetA = numel(setA);
        ar_a = zeros(numsetA,1);
        ac_a = ar_a;
        
        if ~exist('setB','var')
            setB = N2+[1:N*lengthu*N1];
        end
        numsetB = numel(setB);
        
        bi_b = zeros(numsetB,1);
        ar_b = zeros(numsetB,1);
        ac_b = zeros(numsetB,1);
        
        if ~exist('setAA','var') %% adjusted
            setAA = setB(end) + setA;
        end
        numsetAA = numel(setAA);
        ar_aa = zeros(numsetAA,1);
        ac_aa = ar_aa;
        
        %% define the parameters
        
        count = 0;
        count1 = numsetA;
        count2 = numsetA+numsetB; % adjusted
        lenreading = size(d,2)-1;
        % lenreadpone = lenreading+1;
        % lengthu = size(u,1);
        Lenvec = numsetA+numsetB+numsetAA; %Lenvec = N*(N + lengthu*(N1));
        sz = lenreading*N;
        szvec = sz+Lenvec;
        % zn = zeros(1,N); %new
        % znt = zn';
        vecparams = zeros(Lenvec,1);
        meanprior = vecparams;
        % %% Noisy missing preamble
        % diffdata_old = [];
        
        
        %% initialize A, B, vecparams
        % Ap = -rand(N);
        % A = [0 zn ; znt Ap];
        
        %% define A and B and obtain vecparams for reduced system
        if ~exist('A','var')
            A = zeros(N1);
        end
        for i = 1:numsetA
            ar_a(i) = 1+ceil(setA(i)/N);
            ac_a(i) = 1+ N*(1 - ceil(mod(setA(i),N)/N)) + mod(setA(i),N); %ind2sub is slower
            vecparams(i)=A(ar_a(i),ac_a(i));
        end
        
        if ~exist('B','var')
            B = zeros(N1,N1, lengthu);
        end
        
        for i = 1:numsetB
            location = setB(i) -N2;
            bi_b(i) = ceil(location/(N2+N));
            
            check = mod(location,N2+N);
            if check~=0
                ac_b(i) = N1 * (1 - ceil(rem(check,N1)/N1))+rem(check,N1);
                ar_b(i) = 1+ceil(check/N1);
            else
                ar_b(i) = N1;
                ac_b(i) = N1;
            end
            vecparams(i+numsetA) = B(ar_b(i),ac_b(i),bi_b(i));
        end
        
        
        if ~exist('AA','var')
            AA = zeros(N1);
        end
        for i = 1:numsetAA
            location = setAA(i) - N2-N*lengthu*N1;
            ar_aa(i) = 1+ceil(location/N);
            ac_aa(i) = 1+ N*(1 - ceil(mod(location,N)/N)) + mod(location,N); %ind2sub is slower
            vecparams(i+numsetA+numsetB) = AA(ar_aa(i),ac_aa(i));
        end
        %%
        
        %% initialize oldsum
        % oldsum = zeros(N1,N1,lenreading); %
        % Uvec = permute(u,[3 2 1]);
        % parfor j = 1:lenreading %change this to lenreading
        %     oldsum(:,:,j) = sum(bsxfun(@times,Uvec(:,j,:),B),3);
        % end
        % sum_uB = oldsum(:,:,1);
        %% Monitors
        Trend = zeros(Lenvec,30); %?
        Trend(:,1) = vecparams;
        firstorder = zeros(10,1); % new line to monitor first order optimality
        % secderr = 1;
        
        %% Initialize xq
        d = [ones(1,1+ lenreading);d];
        xq = d;
        
        % parfor i = 1:lenreading
        %     store(:,:,i) = expm((pack.t(i+1) -pack.t(i))*(A + oldsum(:,:,i)));
        % end
        
        for i = 1:lenreading
            %     xq(:,i+1) = store(:,:,i)*xq(:,i) ;
            xq(:,i+1) = dopri54c(xq(:,i),delt,A,B,u,i,AA) ;
        end
        
        %% Plot tools %?
        %         figure(1)
        %         plot(pack.t,d(2,:),'b-','Linewidth',2), hold on
        %         hand1 = plot(pack.t,xq(2,:),'k-.','Linewidth',1,'YDataSource','xq(2,:)');
        %         mov1(15)= struct('cdata',[],'colormap',[]);
        %         set(gca,'nextplot','replacechildren')
        %
        %         figure(2)
        %         plot(pack.t,d(3,:),'b-','Linewidth',2), hold on
        %         hand2 = plot(pack.t,xq(3,:),'k-.','Linewidth',1,'YDataSource','xq(3,:)');
        %         mov2 = mov1;
        %         set(gca,'nextplot','replacechildren')
        %
        %         % text position
        %         txx = 0.7*max(pack.t);
        %         txy1 = 0.2*min(d(2,:))+ 0.8*max(d(2,:));
        %         txy2 = 0.2*min(d(3,:))+ 0.8*max(d(3,:));
        %         pos1 = [txx txy1];
        %         pos2 = [txx txy2];
        
        %% M-step variables
        lam = ones(sz,1);
        der = lam;
%         Fstep = zeros(20,1);
        % secder = zeros(sz,sz);
        Cthprior = 10*eye(Lenvec); % larger premultipliers imply less confidence in priors
        %% tight bounds
        % if flag
        % % set Cthprior.
        % end
        
        %% Initialize variables for J matrix/DF
        hwhole = zeros(sz,Lenvec);
        h = zeros(N1,lenreading);
        h1 = h;
        h2=h;
        store1 = h;
        % DF = hwhole(:,1);
        
        span = 2;
        ind = -span:span;
        dex = ind +span+1;
        tbsplinez = zeros(N1,2*span+1);
        % xrange = zeros(1,2*span+1);
        %% Loop utilities
        % outer
        % diffdatamissing = 1;
        % tolmissing = 0.2;
        %inner
        delta = 0.01;
        delta1 = delta;
        tol = 1e-8; % set tolerance
        % tol2 = 0.09;
        conv =1;
        % conv2 =1;
        n = 0;
        %% Loop Flags
        diff = ones(30,1)*Inf;
        % exitflag = 0; %?
        % temp = -1; %?
        % Iter = 30;
        % rcondCemin = 1e-5;
        %% Loop
        % master loop here
        
        while conv > tol
            tic
            %     diffdata=[];
            n=n+1;
            n1 = n+1;
            
            for k = 1:numsetA
                %         ar = 1+ceil(k/N);
                %         ac = 1+ N*(1 - ceil(mod(k,N)/N)) + mod(k,N);
                ar = ar_a(k);
                ac = ac_a(k);
                count = count +1;
                Arac = A(ar,ac);
                parfor i = 1:lenreading
                    Adel=A;
                    tbspline = tbsplinez;
                    xrange = Arac + ind*delta;
                    for g = dex
                        Adel(ar,ac) = xrange(g);
                        %                 tbspline(:,g) = expm((pack.t(i+1) - pack.t(i))*(Adel + oldsum(:,:,i)))*d(:,i) ;
                        tbspline(:,g) = dopri54c(d(:,i),delt,Adel,B,u,i,AA);
                    end
                    store1(:,i) = tbspline(:,span+1);
                    h(:,i) = ppval(fnder(spline(xrange, tbspline)), Arac);
                end
                hwhole(:,count) = reshape(h(2:end,:)',[sz,1]);
            end
            count = 0;
            
            for j = 1:numsetB
                ar = ar_b(j);
                ac = ac_b(j);
                bi = bi_b(j);
                count1 = count1 + 1;
                B_ar_ac_bi = B(ar,ac,bi);
                %         sliced = oldsum(ar,ac,:);
                parfor i = 1:lenreading
                    Bdel=B;
                    %             sum_uB = oldsum(:,:,i);
                    tbspline = tbsplinez;
                    xrange = B_ar_ac_bi + ind*delta1;
                    for g = dex
                        Bdel(ar,ac,bi) = xrange(g);
                        %                 sum_uB(ar,ac) = sliced(i) + u(bi,i)*ind(g)*delta1;
                        %                 tbspline(:,g) = expm((pack.t(i+1) - pack.t(i))*(A + sum_uB))*d(:,i) ;
                        tbspline(:,g) = dopri54c(d(:,i),delt,A,Bdel,u,i,AA)
                    end
                    h1(:,i) =  ppval(fnder(spline(xrange, tbspline)), B_ar_ac_bi) ;
                end
                hwhole(:,count1) = reshape(h1(2:end,:)',[sz,1]);
            end
            count1 = numsetA;
            
            for k = 1:numsetAA
                %         ar = 1+ceil(k/N);
                %         ac = 1+ N*(1 - ceil(mod(k,N)/N)) + mod(k,N);
                ar = ar_aa(k);
                ac = ac_aa(k);
                count2 = count2 +1;
                AArac = AA(ar,ac);
                parfor i = 1:lenreading
                    AAdel=AA;
                    tbspline = tbsplinez;
                    xrange = AArac + ind*delta;
                    for g = dex
                        AAdel(ar,ac) = xrange(g);
                        %                 tbspline(:,g) = expm((pack.t(i+1) - pack.t(i))*(Adel + oldsum(:,:,i)))*d(:,i) ;
                        tbspline(:,g) = dopri54c(d(:,i),delt,A,B,u,i,AAdel);
                    end
                    h2(:,i) = ppval(fnder(spline(xrange, tbspline)), AArac);
                end
                hwhole(:,count2) = reshape(h2(2:end,:)',[sz,1]);
            end
            count2 = numsetA+numsetB;
            
            %% ybar,Jbar
            DF = reshape((d(2:end,2:end)-store1(2:end,:))',[sz,1]);
            
            ybar = [DF; meanprior - vecparams];
            Jbar = [hwhole; eye(Lenvec)];
            
            %% M step
            Cebar = [diag(lam) zeros(sz,Lenvec); zeros(Lenvec,sz)  Cthprior]; % define Cthprior
            P = inv(Cebar) - (Cebar\Jbar)*((Jbar'*(Cebar\Jbar))\Jbar')/Cebar;
            %% verify M step
%             Fstep(n)  = 0.5*(log(abs(det(inv(Cebar)))) - ybar'*(Cebar\ybar)...
%                 - trace(Cthy*Jbar'*(Cebar\Jbar)) + log(abs(det(Cthy))));
            
            %% Continue M Step
            parfor m = 1:sz
                randvar = zeros(szvec,1);
                randvar(m) = sum(P(m,:).*ybar');
                der(m) = -0.5*P(m,m) + 0.5* ybar'*P'*randvar; % check
                randvar(m) = 0;
            end
            secder = -0.5*P.*P';
            secder= secder(1:sz,1:sz);
            
            secderr = rcond(secder);
            if secderr <1e-13
                updatelam = 0;
            else
                updatelam = secder\der;
            end
            %     conv2 = rms(updatelam);
            lam = lam - 0.6*updatelam;
            %% stability monitors
            
            %% Update parameters
            updateparams = (Jbar'*(Cebar\Jbar))\(Jbar'*(Cebar\ybar));
            oldparams = vecparams; %?
            vecparams = oldparams + updateparams;
            
            %% Update A,B. % new
            for j = 1:numsetA
                ar = ar_a(j);
                ac = ac_a(j);
                A(ar,ac) = vecparams(j);
            end
            
            for j = 1:numsetB
                ar = ar_b(j);
                ac = ac_b(j);
                bi = bi_b(j);
                B(ar,ac,bi) = vecparams(j+numsetA);
            end
            
            for j = 1:numsetAA
                ar = ar_a(j);
                ac = ac_a(j);
                AA(ar,ac) = vecparams(j+numsetA+numsetB);
            end
            
            for i = 1:lenreading
                xq(:,i+1) = dopri54c(xq(:,i),delt,A,B,u,i,AA) ;
            end
            diff(n1) = norm(d(2:end,2:end) - xq(2:end,2:end));
            %% Search for best parameters \convergence
            while diff(n1) > diff(n)||isnan(diff(n1))
                %         disp('halve')
                updateparams = 0.6*updateparams;
                vecparams = oldparams + updateparams;
                
                for j = 1:numsetA
                    ar = ar_a(j);
                    ac = ac_a(j);
                    A(ar,ac) = vecparams(j);
                end
                
                for j = 1:numsetB
                    ar = ar_b(j);
                    ac = ac_b(j);
                    bi = bi_b(j);
                    B(ar,ac,bi) = vecparams(j+numsetA);
                end
                
                for j = 1:numsetAA
                    ar = ar_a(j);
                    ac = ac_a(j);
                    AA(ar,ac) = vecparams(j+numsetA+numsetB);
                end
                
                for i = 1:lenreading
                    
                    xq(:,i+1) = dopri54c(xq(:,i),delt,A,B,u,i,AA) ;
                end
                
                diff(n1) = norm(d(2:end,2:end) - xq(2:end,2:end));
                
            end
            
            Trend(:,n1) = vecparams;
            
            firstorder(n) = norm(hwhole,inf);
            
            %% Check conditions for proper exit
            conv = rms(updateparams);
            %% Movies
            %             figure(1)
            %             refreshdata(hand1,'caller')
            %             S2=sprintf('Iter=%d',n);
            %             txmov1 = text('Position',pos1, 'String',S2);
            %             mov1(n) = getscreen(gcf);
            %             delete(txmov1)
            %             drawnow
            %
            %             figure(2)
            %             refreshdata(hand2,'caller')
            %             S3=sprintf('Iter=%d',n);
            %             txmov2 = text('Position',pos2, 'String',S3);
            %             mov2(n) = getscreen(gcf);
            %             delete(txmov2)
            %             drawnow
            
            runtime(n) = toc;
        end
        packD.d   = xq;
        packD.p1  = A;
        packD.p2  = B;
        packD.p3  = AA;
        packD.v   = vecparams;
        packD. C  = Cebar;
        packD.f   = firstorder;
        packD.t   = Trend;
%         packD.do  = pack.do;
        
    end
end

function [y]=dopri54c(y0,h,A,B,u,k,varargin)
t21=0;
y=y0;
a4=[44/45 -56/15 32/9]';
a5=[19372/6561 -25360/2187 64448/6561 -212/729]';
a6=[9017/3168 -355/33 46732/5247 49/176 -5103/18656]';
a7=[35/384 0 500/1113 125/192 -2187/6784 11/84]';
AA1 = varargin{:};
k1=feval(@funcion,t21,y,A,B,u,k,AA1);

k2=feval(@funcion,t21+h/5,y+h*k1/5,A,B,u,k,AA1);
k3=feval(@funcion,t21+3*h/10,y+h*(3*k1+9*k2)/40,A,B,u,k,AA1);
k4=feval(@funcion,t21+4*h/5,y+h*(a4(1)*k1+a4(2)*k2+a4(3)*k3),A,B,u,k,AA1);
k5=feval(@funcion,t21+8*h/9,y+h*(a5(1)*k1+a5(2)*k2+a5(3)*k3+...
    a5(4)*k4),A,B,u,k,AA1);
k6=feval(@funcion,t21+h,y+h*(a6(1)*k1+a6(2)*k2+a6(3)*k3+a6(4)*k4+...
    a6(5)*k5),A,B,u,k,AA1);
y=y+h*(a7(1)*k1+a7(3)*k3+a7(4)*k4+a7(5)*k5+a7(6)*k6);
    function yy = funcion(~,y,A,B,u,k,AA1)
        Uvec = permute(u,[3 2 1]);
        if nargin<7
            yy = A*y +...
                sum(bsxfun(@times,Uvec(:,k,:),B),3)*y;
        else
            yy = A*y +...
                sum(bsxfun(@times,Uvec(:,k,:),B),3)*y+...
                diag(y)*AA1*y;
        end
        
    end
end
