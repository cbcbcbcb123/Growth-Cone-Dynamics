function [LatticeDot] = Plot_lattice_vertical_stiffness(a1, a2, c, b1, b2, b0)
LatticeDot=[];
k=1;
for i=1:c
    x=a1*rand; %格点横坐标
    y=a2*rand; %格点纵坐标
    E=b1+(b2-b1)*y/a2;
    LatticeDot=[LatticeDot;[k,x,y,10^E*10^b0/(10^E+10^b0),0]];
    k=k+1;
end