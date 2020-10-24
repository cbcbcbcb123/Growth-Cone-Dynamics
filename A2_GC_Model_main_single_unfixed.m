clear
clc

Slope=[];
for AA=1:5
ContractionArry_1=[];
VactinArry_1=[];   
for KK=1:1
    KK
% plot grid
L=10000; %宽 ~10000nm = 10um 
H=10000; %高 ~10000nm = 10um 
M=20000; %ECM molecules number, 保持间距为 70nm; 随机分布；
E0=0; %integrin stiffness
E1=-2; %10^-2 ECM stiffness
E2=2; %10^2 ECM stiffness
[Occupy_Function] = Plot_lattice_vertical_stiffness(L, H, M, E1, E2, E0);

% Initiate F-actin
Clutch_Function=[];%提取出覆盖点
L_filaments=2000; %丝状伪足长度 2000 nm
Theta=pi*0/18; %丝状伪足角度
A_top=[4000,4000]; %丝状伪足顶点角度 朝下流动_1（x = 4000nm; y = 4000nm;）； 朝上流动_-1（x = 6000nm; y = 6000nm;）；
Derta_L=35; %丝状伪足宽度 ~75nm
k=1; %第一条
Direction_flow=1; %1_朝下流动；-1_朝上流动；注意顶点要变化
[Occupy_Function, Clutch_Function] = Ini_filaments(L_filaments, Theta, A_top, Derta_L, k, Occupy_Function, Clutch_Function, Direction_flow);

% simulation parameters
MCstep=10000;
t=0;
kon=0.3;
koff=0.1;
Fb=2;
V0=120;
Fstall=200;
Sumforce=0;

L_max=3000;
L_min=1000;

Vpoly=130;

TimeArry=[];
NengageArry=[];
ContractionArry=[];
VactinArry=[];
L_filamentsArry=[];
Num_clutchArry=[];
% run
for cycle=1:MCstep
    % 进度
    if mod(cycle,MCstep/10)==0
        cycle/MCstep
    end
    % Caculate Force-dependent rate
    Num_clutch=size(Clutch_Function,1);
    for j=1:Num_clutch
        if Clutch_Function(j,3)==1
           Clutch_Function(j,5)=koff*exp(Clutch_Function(j,4)*Clutch_Function(j,2)/Fb); 
        end
    end
    % Caculate time
    Tlink=100*ones(Num_clutch,1);
    for j=1:Num_clutch
        if Clutch_Function(j,3)==0
            Tlink(j,1)=-log(rand)/kon;
        else
            Tlink(j,1)=-log(rand)/Clutch_Function(j,5);
        end
    end
    % Min time
    [tmin,X_min_loc]=min(Tlink(:));
    [X_which_link,X_what_happen]=ind2sub(size(Tlink),X_min_loc);
    % Go on
    t=t+tmin;
    % Things happen
    if Clutch_Function(X_which_link,3)==0
        Clutch_Function(X_which_link,3)=1;
    else
        Clutch_Function(X_which_link,3)=0;
        Clutch_Function(X_which_link,4)=0;
        Clutch_Function(X_which_link,5)=100;
    end
    % Update clutch strectch
    Vactin=V0*(1-Sumforce/Fstall);
    for j=1:Num_clutch
        if Clutch_Function(j,3)==1
            Clutch_Function(j,4)= Clutch_Function(j,4)+Vactin*tmin;
        end
    end
    % Sum strectch and bind number
    Sumforce=0;
    Nengage=0;
    for j=1:Num_clutch
        if Clutch_Function(j,3)==1
            Sumforce=Sumforce+Clutch_Function(j,2)*Clutch_Function(j,4);
            Nengage=Nengage+1;
        end
    end
    % Update L 
    L_filaments=L_filaments+(Vpoly-Vactin)*tmin;
    if L_filaments<=L_min || L_filaments>=L_max
        break
    end
    % Add and delete clutch into ABC'D'
    [Occupy_Function, Clutch_Function] = Ini_filaments_add_delete(L_filaments, Theta, A_top, Derta_L, k, Occupy_Function,Clutch_Function, Direction_flow);
    Num_clutch=size(Clutch_Function,1); 
    
    TimeArry(cycle,1)=t;
%     Num_clutchArry(cycle,1)=Num_clutch;
%     NengageArry(cycle,1)=Nengage;
%     ContractionArry(cycle,1)=Sumforce;
%     VactinArry(cycle,1)=Vactin;
    L_filamentsArry(cycle,1)=L_filaments;
end

figure (1)
plot(TimeArry,L_filamentsArry,'o','MarkerSize',2,'MarkerFaceColor',[0,0.75,0.75],'MarkerEdgeColor','none')
hold on
a=polyfit(TimeArry,L_filamentsArry,1);
xi=0:0.001:TimeArry(size(TimeArry,1));
yi=polyval(a,xi);
plot(xi,yi)
a(1)
a(2)

% figure (2)
% plot(TimeArry,NengageArry)
end
Slope=[Slope;a(1)];
end
Slope




