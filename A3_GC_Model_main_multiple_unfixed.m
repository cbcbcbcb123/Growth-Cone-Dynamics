clear
clc

ContractionArry_1=[];
BILI=[];
S1=[];
S2=[];
for KK=0:0.5:0
    KK
    % plot grid
    L=10000; %�� ~10000nm = 10um
    H=20000; %�� ~20000nm = 20um
    M=20000; %ECM molecules number, ���ּ��Ϊ 70nm; ����ֲ���
    E0=0; %integrin stiffness
    E1=-0.5; %10^-2 ECM stiffness
    E2=-0.5; %10^2 ECM stiffness -0.5 for Fig 7
    [Occupy_Function] = Plot_lattice_vertical_stiffness(L, H, M, E1, E2, E0);
    
    % α�㳤�ȣ�A�㣬���
    A_top=[1000,10000]; R_GC=2000; Derta_L=75; L=3000;
    n=21;%���������
    jiaodu=pi/(n-1);
    % α���ռ��λ��
    N_filopodia=zeros(n,6);%��ţ�ռ��1/��ռ��0���Ƕȣ���������A[����]X,Y
    N_filopodia(:,1)=1:n;
    N_filopodia([1,n],3)=pi/2;
    N_filopodia([1,n],5)=A_top(1);
    N_filopodia(1,6)=A_top(2)-R_GC;
    N_filopodia(n,6)=A_top(2)+R_GC;
    for i=2:fix(n/2)
        N_filopodia(i,3)=pi/2+(i-1)*jiaodu;
        N_filopodia(i,5)=A_top(1)-R_GC*cos(pi/2+(i-1)*jiaodu);
        N_filopodia(i,6)=A_top(2)-R_GC*sin(pi/2+(i-1)*jiaodu);
    end
    for i=fix(n/2)+1:n-1
        N_filopodia(i,3)=(i-fix(n/2)-1)*jiaodu;
        N_filopodia(i,5)=A_top(1)+R_GC*cos((i-fix(n/2)-1)*jiaodu);
        N_filopodia(i,6)=A_top(2)-R_GC*sin(pi/2+(i-1)*jiaodu);
    end
    N_filopodia(1:fix(n/2),4)=-1*ones(fix(n/2),1);%1_����������-1_��������
    N_filopodia(fix(n/2)+1:n,4)=ones(n-fix(n/2),1);
    
    %α����Ŀ
    N=21;Count=0;k=[]; %5 for Fig 7
    while Count<N
        R=randi(n);
        if N_filopodia(R,2)==0
            N_filopodia(R,2)=1;
            Count=Count+1;
            k=[k;R];
        end
    end
    
    % Initiate F-actin
    Clutch_Function=[];%��ȡ�����ǵ�
    for i=1:N
        R_top=[N_filopodia(k(i),5),N_filopodia(k(i),6)];
        [Occupy_Function, Clutch_Function] = Ini_filaments(L, N_filopodia(k(i),3), ...
            R_top, Derta_L, N_filopodia(k(i),1),Occupy_Function, Clutch_Function, N_filopodia(k(i),4));
    end
    
%     scatter(Clutch_Function(:,6),Clutch_Function(:,7),10,'filled','MarkerFaceColor','r')
    
    % simulation parameters
    MCstep=5000;
    t=0;
    kon=0.3; %0.15 for Fig 7
    koff=0.1;
    Fb=2;
    V0=120;
    Vpoly1=110;
    Vpoly2=0;
    Fstall=200;%unfixed ���ʼL������ ,3000��Ӧ���ÿ����50��clutch
    Sumforce=zeros(n,1);
    L_filaments=zeros(n,3);%1-���ȣ�2-�仯�ĳ���;3-0/δ��ñ-1��ñ
    L_filaments(k(:),1)=L;
    L_max=6000;
    L_min=10;
    
    Add_pro=0.005; %����filopodia�ĸ���
    Add_cap=0.001; %��ñ�ĸ���
    
    TimeArry=[];
    NengageArry=[];
    ContractionArry=[];
    VactinArry=[];
    L_filamentsArry=[];
    Num_clutchArry=[];
    k_number=[];
    % run
    for cycle=1:MCstep
        % ����
%         if mod(cycle,MCstep/10)==0
%             cycle/MCstep
%         end
        % Caculate Force-dependent rate
        Num_clutch=size(Clutch_Function,1);
        for i=1:Num_clutch
            if Clutch_Function(i,3)==1
                Clutch_Function(i,5)=koff*exp(Clutch_Function(i,4)*Clutch_Function(i,2)/Fb);
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
        Vactin=zeros(n,1);
        for i=1:N
            Vactin(k(i),1)=V0*(1-Sumforce(k(i),1)/Fstall);
        end
        
        for j=1:Num_clutch
            if Clutch_Function(j,3)==1
                Clutch_Function(j,4)= Clutch_Function(j,4)+Vactin(Clutch_Function(j,8),1)*tmin;
            end
        end
        % Sum strectch and bind number
        Sumforce=zeros(n,1);
        Nengage=zeros(n,1);
        
        for j=1:Num_clutch
            if Clutch_Function(j,3)==1
                Sumforce(Clutch_Function(j,8),1)=Sumforce(Clutch_Function(j,8),1)...
                    +Clutch_Function(j,2)*Clutch_Function(j,4);
                Nengage(Clutch_Function(j,8),1)=Nengage(Clutch_Function(j,8),1)+1;
            end
        end
        % Update L
        k_delete=[];
        k2_delete=[];
        for i=1:N
            if L_filaments(k(i),3)==0
                Vpoly=Vpoly1;
            else
                Vpoly=Vpoly2;
            end
%             L_filaments(k(i),1)=L_filaments(k(i),1)+(Vpoly-Vactin(k(i),1))*tmin;
%             L_filaments(k(i),2)=L_filaments(k(i),2)+(Vpoly-Vactin(k(i),1))*tmin;
            L_filaments(k(i),1)=L_filaments(k(i),1)+(Vpoly-Vactin(k(i),1))*0.1;%only for Fig 7
            L_filaments(k(i),2)=L_filaments(k(i),2)+(Vpoly-Vactin(k(i),1))*0.1;
            if abs(L_filaments(k(i),2))>=100 % add and delete clutch into ABC'D'
                R_top=[N_filopodia(k(i),5),N_filopodia(k(i),6)];
                [Occupy_Function, Clutch_Function] = Ini_filaments_add_delete(L_filaments(k(i),1), N_filopodia(k(i),3), ...
                    R_top, Derta_L, N_filopodia(k(i),1),Occupy_Function, Clutch_Function, N_filopodia(k(i),4));
                L_filaments(k(i),2)=0;
            end
            
            if L_filaments(k(i),1)<L_min
                k_delete=[k_delete;k(i)];
                k2_delete=[k2_delete;i];
            end
            
            if max(L_filaments(k(i),1))>L_max
                L_filaments(k(i),1)=L_max;
            end
        end
        
        if size(k_delete,1)>=1
            L_filaments(k_delete,1:3)=zeros(size(k_delete,1),3);
            N_filopodia(k_delete,2)=zeros(size(k_delete,1),1);
        end
        if size(k2_delete,1)>=1
            k(k2_delete,:)=[];
        end
        N=size(k,1);
        
        % Add filopodia
        if N<=n-1 && rand<=Add_pro %�յ�Խ��Խ�������
            M=1;Count=0;
            while Count<M
                R2=randi(n);
                if N_filopodia(R2,2)==0
                    N_filopodia(R2,2)=1;
                    Count=Count+1;
                    k=[k;R2];
                    L_filaments(R2,1)=L;
                end
            end
            
            for i=1:N
                R_top=[N_filopodia(R2,5),N_filopodia(R2,6)];
                [Occupy_Function, Clutch_Function] = Ini_filaments(L, N_filopodia(R2,3), ...
                    R_top, Derta_L, N_filopodia(R2,1),Occupy_Function, Clutch_Function, N_filopodia(R2,4));
            end
        end
        
        % Add cap
        if size(k,1)>=1
            R3=randi(size(k,1));
            Pcap=Add_cap*(1+L_filaments(k(R3),1)/L);
            if rand<=Pcap %���Ǽ���Խ��Խ���׼�ñ��
                if L_filaments(k(R3),3)==0
                    L_filaments(k(R3),3)=1;
                end
            end
        end
        
        TimeArry(cycle,1)=t;
        % Num_clutchArry(cycle,1)=Num_clutch;
        % NengageArry(cycle,:)=Nengage(:,1);
        ContractionArry(cycle,:)=Sumforce(:,1);
%         VactinArry(cycle,:)=Vactin(:,1);
        L_filamentsArry(cycle,:)=L_filaments(:,1);
        k_number(cycle,1)=N;
    end
    
    %�ܽ���һ��
    % ���û�м�ñ��ӻ�����ƽ��״̬������ͼ
    figure (1)
    for i=1:21
        plot(TimeArry,L_filamentsArry(:,i))
        hold on
    end
    % ������̬���Ⱥ�������,��ͼ
%     L_filamentsArry(size(L_filamentsArry,1),:)'
%     mean_ContractionArry=zeros(n,1);
%     for i=1:n
%         mean_ContractionArry(i,1)=mean(ContractionArry(:,i));
%     end
%     mean_ContractionArry
    
    %�ܽ�������
    figure (2)
    subplot(3,1,1)
    plot(TimeArry,k_number)
    
    subplot(3,1,2)
    ContractionArry_0=zeros(MCstep,1);
    for i=1:MCstep
        ContractionArry_0(i,1)=sum(ContractionArry(i,:));
    end
    mean(ContractionArry_0);
    std(ContractionArry_0);
    plot(TimeArry,ContractionArry_0)
    
%     subplot(3,1,3)
%     VactinArry_0=zeros(MCstep,1);
%     for i=1:MCstep
%         VactinArry_0(i,1)=mean(VactinArry(i,:));
%     end
%     plot(TimeArry,VactinArry_0)

    %�ܽ�������
%     figure (3)
%     mean_L_filamentsArry=zeros(n,1);
%     for i=1:n
%         mean_L_filamentsArry(i,1)=mean(L_filamentsArry(:,i));
%     end
%     mean_L_filamentsArry
%     mean(mean_L_filamentsArry(1:10));
%     mean(mean_L_filamentsArry(12:21));
%     BILI=[BILI,mean(mean_L_filamentsArry(1:10))/mean(mean_L_filamentsArry(12:21))]
%     plot(1:n,mean_L_filamentsArry(:,1))
    
%     AA=[];% �ϱ�
%     BB=[];% �±�
%     for i=1:size(Clutch_Function,1)
%         if Clutch_Function(i,7)>A_top(2)
%             AA=[AA;Clutch_Function(i,2)];
%         else
%             BB=[BB;Clutch_Function(i,2)];
%         end
%     end
%     S1=[S1;mean(AA(:))]
%     S2=[S2;mean(BB(:))]
end
