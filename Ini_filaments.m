function [Occupy_Function, Clutch_Function] = Ini_filaments(L_filaments, Theta, A_top, Derta_L, k, Occupy_Function, Clutch_Function, Direction)
%依次输入F-actin长度，角度，一端坐标
%求B点
B_top=[A_top(1,1)+L_filaments*cos(Theta)*Direction,A_top(1,2)+L_filaments*sin(Theta)*Direction];
%求A1，A2，B1，B2
A1=[A_top(1,1)+Derta_L*cos(pi/2+Theta),A_top(1,2)+Derta_L*sin(pi/2+Theta)];
A2=[A_top(1,1)+Derta_L*cos(Theta-pi/2),A_top(1,2)+Derta_L*sin(Theta-pi/2)];
B1=[B_top(1,1)+Derta_L*cos(pi/2+Theta),B_top(1,2)+Derta_L*sin(pi/2+Theta)];
B2=[B_top(1,1)+Derta_L*cos(Theta-pi/2),B_top(1,2)+Derta_L*sin(Theta-pi/2)];
%判断点是否在平行四边形内
xv=[A1(1,1),B1(1,1),B2(1,1),A2(1,1)];
yv=[A1(1,2),B1(1,2),B2(1,2),A2(1,2)];
for i=1:size(Occupy_Function,1)
    if Occupy_Function(i,5)==0
        IN = inpolygon(Occupy_Function(i,2),Occupy_Function(i,3),xv,yv);
        if IN==1
            Occupy_Function(i,5)=k;
            Clutch_Function=[Clutch_Function;Occupy_Function(i,1),Occupy_Function(i,4),0,0,100,Occupy_Function(i,2),Occupy_Function(i,3), k];
        end
    end
end
































