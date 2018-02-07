function [A,B,J,sset]=bnbecon_2edit_(idman, ndes, Nworker, Juu, Jud, Wd, Wn, V, G_hat,Gp,Gd,Gac,Gmvforac, cases, ny,nd)

%ndes is the number of solutions you desire, ie. 10 => returns the global
%top 10 solutions
%mv is the number of items to be selected
%cv is the number of items available to be selected. So, cv>mv
%Nworker is the number of workers, ie processors, available to you to solve
%the problem
























Y = [(Gp/Juu*Jud - Gd)*Wd Wn]; % neglecting implementation error We; this is from Kariwala and Cao
YYT = Y*Y';
% Rank = min([ny+nd, size(Gp,2)]);
% Gscale1=rand(cv,mv);
cv = size(Gp,1);
mv = size(Gp,2);
timer=0;
LB = 0.5;
UB = 2;
ac = 1;
h1j=zeros(1,cv);
timer=0;

%h1j is a vector which you use to for your branching strategy
h1j=1:cv;
flag1=true;
flag2=true;

%F is the fixed set of a node, C is the candidate set of a node, both vectors of logical 1s and 0s
F=false(1,cv);
C=true(1,cv);


% user specifies priority candidate variabels
if ~isempty(idman)
    F(1,idman)=true;
    C(1,idman)=false;
end

%ASTB is your initial incumbent, ie perfection
ASTB=0;
BSTB=0;
% check_lonely_cv =sum(V,2)==0;
% if any(check_lonely_cv)
%     C(1,check_lonely_cv)=false; % remove lonely cvs from candidate
% end


check_single_mv =sum(V,1)==1;
if any(check_single_mv)
    [cvindex_single_mv,~]=find(V(:,check_single_mv));
    F(1,cvindex_single_mv)=true; % add the single cvs to the fixed set
    C(1,cvindex_single_mv)=false;% remove the single cvs from the candidate set
end

num1=nchoosek(sum(C),mv);
OriginalSize=num1;
numpast=num1;
ita=0;
A=Inf(ndes,1);
B=Inf(ndes,1);
J=Inf(ndes,1);
sset=zeros(ndes,mv);






































dumbsolve();

    function dumbsolve()
        while ~isempty(F)
            ita=ita+1;
            Abound=max(A);
            Bbound=max(B);
            %Pruning of supernodes
            if ita>1
                L1=ASTB>Abound;
                L2 = BSTB>Bbound;
                discard=find(L1+L2==2);
                if any(discard)
                    F(discard,:)=[];
                    C(discard,:)=[];
                    ASTB(discard,:)=[];
                    BSTB(discard,:)=[];
                    num1(discard,:)=[];
                end
            end
            
            
            %Tracks Progress
            FC=[];
            CC=[];
            ASTBT=[];
            BSTBT = [];
            nchooseC=[];
            CurrentPos=sum(num1)
            Reduction=-(CurrentPos-numpast)
            PercentSolved=(1-CurrentPos/OriginalSize)*100
            [NumberofNodes,~]=size(F)
            numpast=CurrentPos;
            t=0;
            if NumberofNodes>28000
                flagconst=true;
            else
                flagconst=false;
            end
            
% initialnodes = min(Nworker, nchoosek(cv,mv)); % this should be uncommented if working with small sized cvs
            
            while NumberofNodes<Nworker
                [~,idc]=max(num1);
                fd=F(idc,:);
                cu=C(idc,:);
                fu=fd;
                h1k=h1j;
                h1k(~cu)=Inf;
                [~,id]=min(h1k);
                fu(id)=true;
                cu(id)=false;
                cd = cu;
                
                %%%===== call the subset selection function=====%%%
                cu = subsetsel(fu,V,cu,cv);
                %%%===== end of function call===================%%%
                
                F(idc,:)=[];
                C(idc,:)=[];
                num1(idc,:)=[];
                
                if sum(fu)==mv || sum(fu|cu)==mv
                    flagup = true;
                else
                    flagup = false;
                end
                if sum(fd)==mv || sum(fd|cd)==mv
                    flagdown = true;
                else
                    flagdown = false;
                end
                [a1,b1]=feval(@Gc,fu,cu,ASTB(idc,:),BSTB(idc,:),true,flagup,cases,Gp,G_hat,Juu,YYT,mv,ny,nd); % there is an assumption here that terminal nodes wont be reached in this loop
                [a2,b2]=feval(@Gc,fd,cd,ASTB(idc,:),BSTB(idc,:),false,flagdown,cases,Gp,G_hat,Juu,YYT,mv,ny,nd); % ditto as above
                ASTB(idc,:)=[];
                ASTB=[ASTB ; a1 ; a2];
                BSTB=[BSTB ; b1 ; b2];
                F=[F ; fu ; fd];
                C=[C ; cu ; cd];
                num1=[num1 ; nchoosek(sum(cu),(mv-sum(fu))) ; nchoosek(sum(cd),(mv-sum(fd)))];
                [NumberofNodes,~]=size(F);
            end
            iterat=Inf;
            if ~flagconst
                if NumberofNodes<100*Nworker
                    FC=F;
                    F=[];
                    CC=C;
                    C=[];
                    ASTBT=ASTB;
                    ASTB=[];
                    BSTBT=BSTB;
                    BSTB=[];
                    nchooseC=num1;
                    num1=[];
                else
                    shf=randperm(NumberofNodes,100*Nworker);
                    FC=F(shf,:);
                    F(shf,:)=[];
                    CC=C(shf,:);
                    C(shf,:)=[];
                    ASTBT=ASTB(shf,:);
                    ASTB(shf,:)=[];
                    BSTBT=BSTB(shf,:);
                    BSTB(shf,:)=[];
                    nchooseC=num1(shf,:);
                    num1(shf,:)=[];
                end
            else
                [num1,shf]=sort(num1);
                F=F(shf,:);
                C=C(shf,:);
                ASTB=ASTB(shf,:);
                BSTB=BSTB(shf,:);
                FC=F(1:100*Nworker,:);
                F(1:100*Nworker,:)=[];
                CC=C(1:100*Nworker,:);
                C(1:100*Nworker,:)=[];
                ASTBT=ASTB(1:100*Nworker,:);
                ASTB(1:100*Nworker,:)=[];
                BSTBT=BSTB(1:100*Nworker,:);
                BSTB(1:100*Nworker,:)=[];
                nchooseC=num1(1:100*Nworker,:);
                num1(1:100*Nworker,:)=[];
            end
            [NodestoWorker,~]=size(FC)
            Li=randperm(NodestoWorker);
            FC=FC(Li,:);
            CC=CC(Li,:);
            ASTBT=ASTBT(Li,:);
            BSTBT=BSTBT(Li,:);
            nchooseC=nchooseC(Li,:);
            xx=floor(NodestoWorker/Nworker);
            for n=1:Nworker
                if n==Nworker
                    problems(n).F=FC((n-1)*xx+1:NodestoWorker,:);
                    problems(n).C=CC((n-1)*xx+1:NodestoWorker,:);
                    problems(n).A=ASTBT((n-1)*xx+1:NodestoWorker,:);
                    problems(n).B=BSTBT((n-1)*xx+1:NodestoWorker,:);
                    problems(n).N=nchooseC((n-1)*xx+1:NodestoWorker,:);
                else
                    problems(n).F=FC((n-1)*xx+1:n*xx,:);
                    problems(n).C=CC((n-1)*xx+1:n*xx,:);
                    problems(n).A=ASTBT((n-1)*xx+1:n*xx,:);
                    problems(n).B=BSTBT((n-1)*xx+1:n*xx,:);
                    problems(n).N=nchooseC((n-1)*xx+1:n*xx,:);
                end
            end
            
            %This section is solved in parallel
            parfor z=1:Nworker
                FT=problems(z).F;
                CT=problems(z).C;
                AT=Inf;
                BT=Inf;
                JT=Inf;
                AST=problems(z).A;
                BST=problems(z).B;
                DUMBA=A;
                DUMBB=B;
                DUMBJoint=J;
                ssetT=zeros(1,mv);
                nchooseT=problems(z).N;
                timer=0;
                termin=true;
                toccers=0;
                while termin
                    tic
                    %chooses which node to branch
                    timer=timer+1;
                    if flagconst
                        [~,idn]=min(nchooseT);
                    else
                        if toccers<90
                            [~,idn]=min(AST.*BST);
                        else
                            [~,idn]=min(nchooseT);
                        end
                    end
                    fd=FT(idn,:);
                    cu=CT(idn,:);
                    cd = cu;
                    AS=AST(idn,:);
                    BS=BST(idn,:);
                    fu=fd;
                    
                    %chooses which variable to branch if it isnt a terminal
                    %node already
                    if sum(fu)<mv
                        h1k=h1j;
                        h1k(~cu)=Inf;
                        [~,id]=min(h1k);
                        fu(id)=true;
                        cu(id)=false;
                    end
                    
                    if sum(fd|cd)>mv
                        cd=cu;
                    end
                    
                    %removing node to be branched
                    FT(idn,:)=[];
                    CT(idn,:)=[];
                    nchooseT(idn,:)=[];
                    AST(idn,:)=[];
                    BST(idn,:)=[];
                    flag1=true;
                    flag2=true;
                    
                    %%%===== call the subset selection function=====%%%
                    try
                    cu = subsetsel(fu,V,cu,cv);
                    catch ME
                        if strcmp(ME.identifier,'MATLAB:badsubscript')
                            cu
                            fu
                            idn
                        end
                    end
                        
                    %%%===== end of function call===================%%%
                    
                    %checking whether node1 is terminal
                    cc=sum(cu);
                    k=mv-sum(fu);
                    if cc<k
                        flag1=false % this should never be reached
                    end
                    
                    if sum(cd)<mv-sum(fd)
                        flag2 = false % this should never be reached
                    end
                    
                    if (k==0||k==cc)&&flag1
                        [a,b]=feval(@Gc,fu,cu,0,0,true,true,cases,Gp,G_hat,Juu,YYT,mv,ny,nd);
                        Abound=min([max(AT) max(DUMBA)]);
                        Bbound=min([max(BT) max(DUMBB)]);
                        Jointbound = min([max(JT) max(DUMBJoint)]);
                        jscale = a*b;
                        
                        %%%===== call the RGA function=====%%%
                        %  RGAd_ac = getRGA(k,fu,cu,Gp,Gac,Gmvforac,ac); % this obtains the RGA of diagonal elements for active constraints
                        %%%===== end of function call======%%%
                        RGAd_ac = 1; % used for the test cases where i wasnt handling rga
                        
%                         try 
%                            (jscale<=Jointbound &&(any(RGAd_ac>=LB) && any(RGAd_ac<=UB)))
%                         catch 
%                             jscale
%                             Jointbound
%                         end
                        if (jscale<=Jointbound &&(any(RGAd_ac>=LB) && any(RGAd_ac<=UB)))
                            AT=[a ; AT];
                            BT=[b ; BT];
                            DUMBA = [a ; DUMBA];
                            DUMBB = [b ; DUMBB];
                            JT = [jscale ; JT];
                            if k==0
                                ssetT=[find(fu) ; ssetT];
                            else
                                s=(fu|cu);
                                ssetT=[find(s) ; ssetT];
                            end
                            if size(AT,1)>ndes
                                [~,IDDD]=max(JT);
                                AT(IDDD,:)=[];
                                BT(IDDD,:) = [];
                                JT(IDDD,:) = [];
                                ssetT(IDDD,:)=[];
                            end
                        end
                        flag1=false;
                    end
                    
                    %checking whether node2 is terminal
                    cc=sum(cd);
                    k=mv-sum(fd);
%                     try
%                        (k==0||k==cc)&&flag2
%                     catch ME
%                         if strcmp(ME.identifier,'MATLAB:nonLogicalConditional')
%                             cd
%                             fd
%                             fu
%                             cu
%                             idn
%                             AST
%                             BST
%                         end
%                     end
                    if (k==0||k==cc)&&flag2
                        [a,b]=feval(@Gc,fd,cd,0,0,true,true,cases,Gp,G_hat,Juu,YYT,mv,ny,nd);
                        Abound=min([max(AT) max(DUMBA)]);
                        Bbound=min([max(BT) max(DUMBB)]);
                        Jointbound = min([max(JT) max(DUMBJoint)]);
                        jscale = a*b;
                        
                                                %%%===== call the RGA function=====%%%
%                         RGAd_ac = getRGA(k,fd,cd,Gp,Gac,Gmvforac,ac); % this obtains the RGA of diagonal elements for active constraints
                        %%%===== end of function call======%%%
                        RGAd_ac = 1;% used for the test cases where i wasnt handling rga
                        
                        if  (jscale<=Jointbound &&(any(RGAd_ac>=LB) && any(RGAd_ac<=UB)))
                            AT=[a ; AT];
                            BT=[b ; BT];
                            DUMBA = [a ; DUMBA];
                            DUMBB = [b ; DUMBB];
                            JT = [jscale ; JT];
                            if k==0
                                ssetT=[find(fd) ; ssetT];
                            else
                                s=(fd|cd);
                                ssetT=[find(s) ; ssetT];
                            end
                            if size(AT,1)>ndes
                                [~,IDDD]=max(JT);
                                AT(IDDD,:)=[];
                                BT(IDDD,:)=[];
                                JT(IDDD,:) = [];
                                ssetT(IDDD,:)=[];
                            end
                        end
                        flag2=false;
                    end
                    
                    %Pruning of Branch 1
                    if flag1
                        [a,b]=feval(@Gc,fu,cu,AS,BS,true,false,cases,Gp,G_hat,Juu,YYT,mv,ny,nd);
                        Abound=min([max(AT) max(DUMBA)]);
                        Bbound=min([max(BT) max(DUMBB)]);
                        Jointbound = min([max(JT) max(DUMBJoint)]);
                        jscale = a*b;
                        
                        if a<=Abound || b<=Bbound||jscale<=Jointbound
                            AST=[AST ; a];
                            BST=[BST ; b];
                            FT=[FT ; fu];
                            CT=[CT ; cu];
                            nchooseT=[nchooseT ; nchoosek(sum(cu),mv-sum(fu))];
                        else
                            flag1=false;
                        end
                    end
                    
                    %Pruning of Branch 2
                    if flag2
                        [a,b]=feval(@Gc,fd,cd,AS,BS,false,false,cases,Gp,G_hat,Juu,YYT,mv,ny,nd);
                        Abound=min([max(AT) max(DUMBA)]);
                        Bbound=min([max(BT) max(DUMBB)]);
                        Jointbound = min([max(JT) max(DUMBJoint)]);
                        jscale = a*b;
                        
                        if a<=Abound || b<=Bbound||jscale<=Jointbound
                            AST=[AST ; a];
                            BST=[BST ; b];
                            FT=[FT ; fd];
                            CT=[CT ; cd];
                            nchooseT=[nchooseT ; nchoosek(sum(cd),mv-sum(fd))];
                        else
                            flag2=false;
                        end
                    end
                    dtt=toc;
                    toccers=toccers+dtt;
                    if timer>iterat || isempty(FT) || toccers>300
                        termin=false;
                    end
                end
                
                %Storing existing nodes and solutions from iteration
                [ttte,~]=size(FT);
                if ttte~=0
                    F=[F ; FT];
                    C=[C ; CT];
                    ASTB=[ASTB ; AST];
                    BSTB=[BSTB ; BST];
                    num1=[num1;nchooseT];
                end
                SOLUTIONS(z).D1=ssetT;
                SOLUTIONS(z).D2=AT;
                SOLUTIONS(z).D3=BT;
                SOLUTIONS(z).D4=JT;
            end
            
            %Post-processing of solutions
            for i=1:Nworker
                sset=[sset ; SOLUTIONS(i).D1];
                A=[A ; SOLUTIONS(i).D2];
                B=[B ;  SOLUTIONS(i).D3];
                J=[J ;  SOLUTIONS(i).D4];
            end
            
            if size(sset,1)>ndes
                [J,idpd]=sort(J);
                A=A(idpd,:);
                B=B(idpd,:);
                sset=sset(idpd,:);
                idsss=size(A);
                A((ndes+1):idsss,:)=[];
                B((ndes+1):idsss,:)=[];
                J((ndes+1):idsss,:)=[];
                sset((ndes+1):idsss,:)=[];
            end
            SOLUTIONS=[];
        end
    end
end
%% functions
% obtain RGA
function RGAd_ac = getRGA(k,f,c,Gp,Gac,Gmvforac,ac)

if k==0
    Gpry = Gp(f,:);
    Gm = Gmvforac(f);
else
    Gpry = Gp(f|c,:);
    Gm = Gmvforac(f|c);
end
Gall = [Gac;Gm Gpry]; % Gac is the gain matrix for the active constraint, Gmvfor ac are the corresponding columns Gpry is the subset of the gain matrix chosen
RGA = Gall.*inv(Gall');
RGAd = diag(RGA);
RGAd_ac = RGAd(1:ac);
end
% subset selection
function c2 = subsetsel(f2,V,c2,cv)
temp = logical(f2);
Vt=V(temp,:);
Vtt=V;
temp=sum(Vt,2);
L=temp==1; % check if some cvs exist that are only controlled by one mv

[~,temp6] = find(Vt(L,:)); % obtain the index of the manipulated variables which cvs with only one mvs are being controlled by
temp7=1:size(Vt,2);
temp7(:,temp6)=[];
temp10=sum(Vtt(:,temp6),2);
temp11=false(cv,1);
temp11(logical(temp10),1)=true;

temp11=temp11-sum(Vtt(:,temp7),2); % check other cvs being controlled by the mvs responsible for the control of cvs with only one mv
S=temp11==1;
if any(S)
    c2(1,S)=false; % remove these ones from the candidates
end
end


%This is your objective function
%flag2 defines rather the node is terminal or not, logical 1 or 0
%flag defines rather the node has been upwardly branched or downwardly branched, 1 == upward branched, 0 == downward branched
function [a,b]=Gc(f,c,as,bs,flag,flag2,cases,Gp,G_hat,Juu,YYT,mv,ny,nd)
f = logical(f);
c = logical(c);
k=mv-sum(f);
warning('off', 'all')
if flag2
    if k==0
        a = upordown(f);
        b = min(nonzeros(svd(G_hat(f,:))))^-1;
    else
        a = upordown((f|c));
        b = (min(nonzeros(svd(G_hat(f|c,:)))))^-1;
    end
else
    if flag
        a=max(upordown(f),as);
        b = max(min(nonzeros(svd(G_hat(f,:))))^-1,bs); 
    else   
        a= max(upordown((f|c)),as);
        b = max((min(nonzeros(svd(G_hat(f|c,:)))))^-1,bs);
    end
end
    function a = upordown(x)
        x = logical(x);
        Gtilde = Gp/(chol(Juu));
        
%         try
%             R = chol(YYT(x,x));
%             fun = eig(R'\Gtilde(x,:)*Gtilde(x,:)'/R);
%         catch ME
%             if strcmp(ME.identifier,'MATLAB:badsubscript')
%                 x = logical(x);
%                 R = chol(YYT(x,x));
%                 fun = eig(R'\Gtilde(x,:)*Gtilde(x,:)'/R);
%             elseif strcmp(ME.identifier,'MATLAB:posdef')
%                 fun = eig(Gtilde(x,:)'*(YYT(x,x)\Gtilde(x,:)));
%             else
%                 ME
%             end
%         end
        if sum(x)<mv
            R = chol(YYT(x,x));
            fun = svd(R'\Gtilde(x,:)*Gtilde(x,:)'/R);
        else
            fun = svd(Gtilde(x,:)'*(YYT(x,x)\Gtilde(x,:)));
        end

        
        if sum(x)<=mv % upwardly branched
            fun = fun(1:sum(x));
        else % downwardly branched
            fun = fun(1:mv);
        end
        
%         if any(fun<0)
%             disp('negative eig:') % this should never happen
%             fun(fun<0)
%             find(x)
%         end
        
        if strcmp(cases,'worst')
            a = 0.5/(min(fun));
        elseif strcmp(cases,'average')
%             a = 1/6/(ny+nd)*sqrt(sum(1./fun));
            a = (1/6/mv)*(sum(1./fun));  % redefined the average loss to my derivation.
        else
            error('undefined case: use "worst" or "average"')
        end
    end

end
