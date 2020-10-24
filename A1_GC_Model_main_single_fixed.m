clear
clc

Ave_stiffness=[];
for aa=-2:0.5:2
    aa
ContractionArry_1=[];
VactinArry_1=[];
Clutch_Function_1=[];
Num_clutchArry_1=[];
Clutch_FunctionArry_1=[];
for KK=1:10
% plot grid
L=10000; %宽 ~10000nm = 10um 
H=10000; %高 ~10000nm = 10um 
M=20000; %ECM molecules number, 保持间距为 70nm; 随机分布；
E0=0; %integrin stiffness
E1=-2; %10^-2 ECM stiffness
E2=aa; %10^2 ECM stiffness
[Occupy_Function] = Plot_lattice_vertical_stiffness(L, H, M, E1, E2, E0);
% figure (2)
% scatter3(Occupy_Function(:,2),Occupy_Function(:,3),Occupy_Function(:,4))

% Initiate F-actin
Clutch_Function=[];%提取出覆盖点
L_filaments=8000; %丝状伪足长度 1000~8000 nm
Theta=pi*0/18; %丝状伪足角度
A_top=[500,500]; %丝状伪足顶点角度 x = 500nm; y = 500nm;
Derta_L=35; %丝状伪足宽度 ~75nm
k=1; %第一条
Direction_flow=1; %1_朝下流动；-1_朝上流动
[Occupy_Function, Clutch_Function] = Ini_filaments(L_filaments, Theta, A_top, Derta_L, k, Occupy_Function, Clutch_Function, Direction_flow);

% Num_clutch=size(Clutch_Function,1);
% hold on
% for i=1:Num_clutch
% scatter(Occupy_Function(Clutch_Function(i,1),2),...
%     Occupy_Function(Clutch_Function(i,1),3),10,'filled','MarkerFaceColor','r')
% end

% simulation parameters
MCstep=100000;
t=0;
kon=0.3;
koff=0.1;
Fb=2;
V0=120;
Fstall=200;
Sumforce=0;

TimeArry=[];
ContractionArry=[];
VactinArry=[];
Num_clutchArry=[];
Clutch_FunctionArry=[];
% run
for cycle=1:MCstep
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
    for j=1:Num_clutch
        if Clutch_Function(j,3)==1
            Sumforce=Sumforce+Clutch_Function(j,2)*Clutch_Function(j,4);
        end
    end

    TimeArry(cycle,1)=t;
    Num_clutchArry(cycle,1)=Num_clutch;
    ContractionArry(cycle,1)=Sumforce;
    Clutch_FunctionArry(cycle,1)=mean(Clutch_Function(:,2));
%     VactinArry(cycle,1)=Vactin;
end

% figure (1)
% 判断是否平衡
% plot(TimeArry,ContractionArry)

% 计算平均值
ContractionArry_1=[ContractionArry_1;mean(ContractionArry)];
Num_clutchArry_1=[Num_clutchArry_1;mean(Num_clutchArry)];
Clutch_FunctionArry_1=[Clutch_FunctionArry_1;mean(Clutch_FunctionArry)];
% VactinArry_1=[VactinArry_1;mean(VactinArry)];

% 输出有效刚度值
% format longG
% Clutch_Function(:,2)
% Clutch_Function_1=[Clutch_Function_1;mean(Clutch_Function(:,2))];
end
format longG
% ContractionArry_1
% Num_clutchArry_1
% Clutch_FunctionArry_1
Ave_stiffness=[Ave_stiffness;mean(ContractionArry_1)];
end
Ave_stiffness





