%% This MAtlab script estimates parameters of A,B,C matrices from the following bilinear model
% \dot{x} = Ax + \sum_{j=1}B_ju_jx + Cu
%using the expectation maximization algorithm.
% data must be provided where size(data) = [N, lenreading]
% u must also be provided where size(u) = [lengthu, lenreading]
% if some elements of matrix A are disallowed from taking non zero values, the set of indices must be provided in setA 
% same follows for setB

%% ismember check
% idx = arrayfun(@(x)find(A==x,1),B);
runtime = zeros(20,1);
%% define the set of nonzero parameters in set A/B
lengthu = size(u,1);
N = size(data,1);
N1 = N+1;
N2 = N^2;

if ~exist('setA','var')
    setA = 1:N2;
    flag = true;
end

numsetA = numel(setA);
ar_a = zeros(numsetA,1);
ac_a = ar_a;

if ~exist('setB','var') 
    setB = N2+[1:N*size(u,1)*N1];
end 
numsetB = numel(setB);

bi_b = zeros(numsetB,1);
ar_b = zeros(numsetB,1);
ac_b = zeros(numsetB,1);

%% define the parameters

count = 0;
count1 = numsetA;
lenreading = size(data,2)-1;
lenreadpone = lenreading+1;
% lengthu = size(u,1);
Lenvec = numsetA+numsetB; %Lenvec = N*(N + lengthu*(N1));
sz = lenreading*N;
szvec = sz+Lenvec;
zn = zeros(1,N); %new
znt = zn';
vecparams = zeros(Lenvec,1);
meanprior = vecparams;
% %% Noisy missing preamble
diffdata_old = [];


%% define A and B and obtain vecparams for reduced system
A = zeros(N1);
for i = 1:numsetA
    ar_a(i) = 1+ceil(setA(i)/N);
    ac_a(i) = 1+ N*(1 - ceil(mod(setA(i),N)/N)) + mod(setA(i),N); %ind2sub is slower
    A(ar_a(i),ac_a(i)) = rand; 
    vecparams(i) = A(ar_a(i),ac_a(i));
end
Adel = A;

B = zeros(N1,N1, lengthu);

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
    B(ar_b(i),ac_b(i),bi_b(i)) = rand;
    vecparams(i+numsetA) = B(ar_b(i),ac_b(i),bi_b(i));
end

%% initialize oldsum
oldsum = zeros(N1,N1,lenreading); %change this to lenreading
Uvec = permute(u,[3 2 1]);
parfor j = 1:lenreading %change this to lenreading
    oldsum(:,:,j) = sum(bsxfun(@times,Uvec(:,j,:),B),3);
end
sum_uB = oldsum(:,:,1);
%% Monitors
Trend = zeros(Lenvec,30); %?
firstorder = zeros(10,1); % new line to monitor first order optimality
secderr = 1;
% secderr = ones(30,1);
% Cthyr = secderr;
% Cer = secderr;

%% Initialize xq
data = [ones(1,1+ lenreading);data];
xq = data;

parfor i = 1:lenreading
    store(:,:,i) = expm((t(i+1) -t(i))*(A + oldsum(:,:,i)));
end

for i = 1:lenreading
    xq(:,i+1) = store(:,:,i)*xq(:,i) ;
end

%% Plot tools %?
% figure(1)
% plot(t,data(2,:),'b-','Linewidth',2), hold on
% hand1 = plot(t,xq(2,:),'k-.','Linewidth',1,'YDataSource','xq(2,:)');
% mov1(15)= struct('cdata',[],'colormap',[]);
% set(gca,'nextplot','replacechildren')
%
% figure(2)
% plot(t,data(3,:),'b-','Linewidth',2), hold on
% hand2 = plot(t,xq(3,:),'k-.','Linewidth',1,'YDataSource','xq(3,:)');
% mov2 = mov1;
% set(gca,'nextplot','replacechildren')
%
% % text position
% txx = 0.7*max(t);
% txy1 = 0.2*min(data(2,:))+ 0.8*max(data(2,:));
% txy2 = 0.2*min(data(3,:))+ 0.8*max(data(3,:));
% pos1 = [txx txy1];
% pos2 = [txx txy2];

%% M-step variables
lam = ones(sz,1);
der = lam;
secder = zeros(sz,sz);
Cthprior = 20*eye(Lenvec); % larger premultipliers imply less confidence in priors
%% tight bounds

if flag
% tight_bounds = [20,21,25,26];
% elem_num = zeros(size(tight_bounds));
% parfor i = 1:numel(tight_bounds)
%     elem_num(i) = tight_bounds(i)*(Lenvec+1)-Lenvec;
% end
% %
% elem_num1 = zeros(size(tight_bounds2));
% parfor i = 1:numel(tight_bounds2)
%     elem_num1(i) = tight_bounds2(i)*(Lenvec+1)-Lenvec;
% end
% elem_num2 = zeros(size(tight_bounds3));
% parfor i = 1:numel(tight_bounds3)
%     elem_num2(i) = tight_bounds3(i)*(Lenvec+1)-Lenvec;
% end
% Cthprior(elem_num) = 0.01;
% Cthprior(elem_num1) = 0.001;
% Cthprior(elem_num2) = 0.01;
end

%% Initialize variables for J matrix
hwhole = zeros(sz,Lenvec);
h = zeros(N1,lenreading);
h1 = zeros(N1,lenreading);
DF = hwhole(:,1);

span = 2;
ind = -span:span;
dex = ind +span+1;
tbsplinez = zeros(N1,2*span+1);
xrange = zeros(1,2*span+1);
%% Loop utilities
% outer
diffdatamissing = 1;
tolmissing = 0.2;
%inner
delta = 0.01;
delta1 = delta;
tol = 1e-8; % set tolerance
tol2 = 0.09;
conv =1;
conv2 =1;
n = 0;
%% Loop Flags
diff = ones(30,1)*Inf;
exitflag = 0; %?
temp = -1; %?
Iter = 30;
% rcondCemin = 1e-5;
%% Loop
% master loop here

while conv > tol
    tic
    diffdata=[];
    n=n+1
    n1 = n+1;
    %         if n<=3
    %             prev = data;
    %         else
    %             prev = xq;
    %         end

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
                tbspline(:,g) = expm((t(i+1) - t(i))*(Adel + oldsum(:,:,i)))*data(:,i) ;
            end
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
        sliced = oldsum(ar,ac,:);
        parfor i = 1:lenreading
            sum_uB = oldsum(:,:,i);
            tbspline = tbsplinez;
            xrange = B_ar_ac_bi + ind*delta1;
            for g = dex
                sum_uB(ar,ac) = sliced(i) + u(bi,i)*ind(g)*delta1;
                tbspline(:,g) = expm((t(i+1) - t(i))*(A + sum_uB))*data(:,i) ;
            end
            h1(:,i) =  ppval(fnder(spline(xrange, tbspline)), B_ar_ac_bi) ;
        end
        hwhole(:,count1) = reshape(h1(2:end,:)',[sz,1]);
    end
    count1 = numsetA;

    
    %% ybar,Jbar
    %         tic
    for i=2:lenreadpone
        df = data(:,i)- store(:,:,i-1)*data(:,i-1); % data row is y, col is time
        DF(i-1:lenreading:(N-1)*lenreading +i-1) = df(2:end); % em for missing data could come here
    end
    %         toc
    
    %         tic
    %         parfor i=1:lenreading
    %             df(:,i) = data(:,i+1)- store(:,:,i)*data(:,i);
    %         end
    %         DF = reshape(df(2:end,:)',[sz 1]);
    %         toc
    
    ybar = [DF; meanprior - vecparams];
    Jbar = [hwhole; eye(Lenvec)];
    
    %% M step
    Cebar = [diag(lam) zeros(sz,Lenvec); zeros(Lenvec,sz)  Cthprior]; % define Cthprior
    P = inv(Cebar) - (Cebar\Jbar)*((Jbar'*(Cebar\Jbar))\Jbar')/Cebar;
    %% verify M step
    %     F(n)  = 0.5*(log(abs(det(inv(Cebar)))) - ybar'*(Cebar\ybar)...
    %         - trace(Cthy*Jbar'*(Cebar\Jbar)) + log(abs(det(Cthy))));
    
    %% Continue M Step
    %         while conv2 > tol2
    
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
        updatelam = 0
        %             break
    else
        updatelam = secder\der;
    end
    conv2 = rms(updatelam);
    lam = lam - 0.6*updatelam;
    %% stability monitors
    
    %% Update parameters
    % updateparams = (Jbar'*(Cebar\Jbar))\(Jbar'*(Cebar\ybar));
    updateparams = (Jbar'*(Cebar\Jbar))\(Jbar'*(Cebar\ybar));
    oldparams = vecparams; %?
    vecparams = oldparams + updateparams;
    
    %% Update A,B. % new
    %     Apre = vec2mat(vecparams(1:N2),N);
    %     A = [0 zn ; znt Apre];
    for j = 1:numsetA
        ar = ar_a(j);
        ac = ac_a(j);
        A(ar,ac) = vecparams(j);
    end
    
    for j = 1:numsetB
        %         jj = j-numsetA;
        ar = ar_b(j);
        ac = ac_b(j);
        bi = bi_b(j);
        B(ar,ac,bi) = vecparams(j+numsetA);
    end
    %     Bpre = vec2mat(vecparams(N2+1:end),N1);
    %
    %     for k = 1:lengthu
    %         B(2:end,:,k) = Bpre(1+(k-1)*N: k*N,:);
    %     end
    
    parfor j = 1:lenreading
        oldsum(:,:,j) = sum(bsxfun(@times,Uvec(:,j,:),B),3);
    end
    
    parfor i = 1:lenreading
        store(:,:,i) = expm((t(i+1) -t(i))*(A + oldsum(:,:,i)));
    end
    
    for i = 1:lenreading
        xq(:,i+1) = store(:,:,i)*xq(:,i) ;
    end
    diff(n1) = norm(data(2:end,2:end) - xq(2:end,2:end))
    %% Search for best parameters \convergence
    while diff(n1) > diff(n)
        disp('halve')
        updateparams = 0.6*updateparams;
        vecparams = oldparams + updateparams;
        
        for j = 1:numsetA
            ar = ar_a(j);
            ac = ac_a(j);
            A(ar,ac) = vecparams(j);
        end
        
        for j = 1:numsetB
            %         jj = j-numsetA;
            ar = ar_b(j);
            ac = ac_b(j);
            bi = bi_b(j);
            B(ar,ac,bi) = vecparams(j+numsetA);
        end
        
        parfor j = 1:lenreading
            oldsum(:,:,j) = sum(bsxfun(@times,Uvec(:,j,:),B),3);
        end
        
        parfor i = 1:lenreading
            store(:,:,i) = expm((t(i+1) -t(i))*(A + oldsum(:,:,i)));
        end
        
        for i = 1:lenreading
            xq(:,i+1) = store(:,:,i)*xq(:,i) ;
        end
        
        diff(n1) = norm(data(2:end,2:end) - xq(2:end,2:end))
        
    end
    Adel = A;
    sum_uB = oldsum(:,:,1);
    Trend(:,n) = vecparams;
    
    firstorder(n) = norm(hwhole,inf);
    %%
    %         for i = 1: N
    %             k = find(ind(i,:)==0,1,'first')-1;
    %             data(i+1,ind(i,1:k)) = xq(i+1,ind(i,1:k));
    %             diffdata = [diffdata; xq(i+1,ind(i,1:k))'];
    %         end
    
    %     if diff(n1)> min(diff) && diff(n)> min(diff)
    %     if n > 3 && diff(n1)> min(diff)
    %         exitflag = -2;% variable unused for all previous iterations but the last
    %         vecparams = Trend(:,n-1);
    %         break
    %     end
    
    %% Check conditions for proper exit
    conv = rms(updateparams);
    %% Movies
    %         figure(1)
    %         refreshdata(hand1,'caller')
    %         S2=sprintf('Iter=%d',n);
    %         txmov1 = text('Position',pos1, 'String',S2);
    %         mov1(n) = getscreen(gcf);
    %         delete(txmov1)
    %         drawnow
    %
    %         figure(2)
    %         refreshdata(hand2,'caller')
    %         S3=sprintf('Iter=%d',n);
    %         txmov2 = text('Position',pos2, 'String',S3);
    %         mov2(n) = getscreen(gcf);
    %         delete(txmov2)
    %         drawnow
    runtime(n) = toc;
end
