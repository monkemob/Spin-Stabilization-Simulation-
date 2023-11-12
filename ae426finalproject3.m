clc
clear
close all
%lvlh to eci
%initial q 
time=[0:50]; %s
Qle=[0 0 1; cosd(30) -sind(30) 0; sind(30) cosd(30) 0];
[V,D]=eig(Qle);
%D(:,1)=1 0 0 therfore V(:,1) is the axis of rotation
rotangle=acos(.189); %rad
uhat=V(:,1);
qinitial=[cos(rotangle/2);sin(rotangle).*uhat];
%lvlh to eci
[T,X]=ode45(@diffeq,time,[0;.00113;0;qinitial(1);qinitial(2);qinitial(3);qinitial(4)]);
%body to lvlh
%ae313 calculations taken from a previously written code
mu=398600;
a=6778;
f=[0:1:51];
e=0;
i=30;
w=0;
omega=0;
p=a*(1-e^2);
if p==0
    p=a;
end
h=sqrt(p*mu);
r=p./(1+e.*cos(f));
rpqw=[r.*cos(f); r.*sin(f); zeros(1,length(f))];
vpqw=(mu/h).*[-sin(f); e+cos(f); zeros(1,length(f))];
v=[];
for k=1:length(f)
v(k)=norm(vpqw(:,k));
end

rot=[cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1]*[1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)]*[cos(w) -sin(w) 0; sin(w) cos(w) 0;0 0 1];
%for reference:
%rot1=[cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1];
%rot2=[1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
%rot3=[cos(w) -sin(w) 0; sin(w) cos(w) 0;0 0 1];

for o=1:length(f)
vvec(:,o)=rot*vpqw(:,o);
rvec(:,o)=rot*rpqw(:,o);
end

for j=1:51
rlvlh(1,j)=rvec(1,j).*(X(j,4)^2-X(j,7)^2-X(j,6)^2+X(j,5)^2)+rvec(2,j).*(2*X(j,5)*X(j,6)-2*X(j,4)*X(j,7))+rvec(3,j).*2*X(j,4)*X(j,6);
rlvlh(2,j)=rvec(1,j).*(2*X(j,6)*X(j,5)-2*X(j,4)*X(j,7))+rvec(2,j).*(X(j,6)^2-X(j,7)^2+X(j,4)^2-X(j,5)^2)+rvec(3,j).*(2*X(j,4)*X(j,5)-2*X(j,6)*X(j,7));
rlvlh(3,j)=rvec(1,j).*2*X(j,4)*X(j,6)+rvec(2,j).*(2*X(j,6)*X(j,7)-2*X(j,4)*X(j,5))+rvec(3,j).*(-X(j,7)^2-X(j,6)^2+X(j,5)^2+X(j,4)^2);
end
rlvlh=[zeros(1,51);rlvlh]; %from def of q multiplication

[T2,X2]=ode45(@diffeq2,time,[0;0;0;1;0;0;0]);

%lvlh to eci (will be using opposite)
function dydt=diffeq(t,y) 
dy1=0;
dy2=0;
dy3=0;
dy4=.5*(-y(5)*y(1)-y(6)*y(2)-y(7)*y(3));
dy5=.5*(y(4)*y(1)-y(7)*y(2)+y(6)*y(3));
dy6=.5*(y(7)*y(1)+y(4)*y(2)-y(5)*y(3));
dy7=.5*(-y(6)*y(1)+y(5)*y(2)+y(4)*y(3));
dydt=[dy1;dy2;dy3;dy4;dy5;dy6;dy7];
end
%for body to lvlh disturbance
function dydt2=diffeq2(t,y) 
%angular velocity
dy1=1+3*mu/(6558*10^3)^5.*(rlvlh(2,1).*(2*y(6)*y(5)-2*y(4)*y(7))+rlvlh(3,1).*(y(6)^2-y(7)^2+y(4)^2-y(5)^2)+rlvlh(4,1).*(2*y(4)*y(5)-2*y(6)*y(7))).*(rlvlh(2,1).*2*y(4)*y(6)+rlvlh(3,1).*(2*y(6)*y(7)-2*y(4)*y(5))+rlvlh(4,1).*(-y(7)^2-y(6)^2+y(5)^2+y(4)^2))-2*y(2)*y(3);
dy2=1+1.5*y(1)*y(3);
dy3=2+3*mu/(6558*10^3)^5.*(rlvlh(2,1).*(2*y(6)*y(5)-2*y(4)*y(7))+rlvlh(3,1).*(y(6)^2-y(7)^2+y(4)^2-y(5)^2)+rlvlh(4,1).*(2*y(4)*y(5)-2*y(6)*y(7))).*(rlvlh(2,1).*(y(4)^2-y(7)^2-y(6)^2+y(5)^2)+rlvlh(3,1).*(2*y(5)*y(6)-2*y(4)*y(7))+rlvlh(4,1).*2*y(4)*y(6))-.25*y(1)*y(2);
%qs. this format is same for all diffeq functions
dy4=.5*(-y(5)*y(1)-y(6)*y(2)-y(7)*y(3));
dy5=.5*(y(4)*y(1)-y(7)*y(2)+y(6)*y(3));
dy6=.5*(y(7)*y(1)+y(4)*y(2)-y(5)*y(3));
dy7=.5*(-y(6)*y(1)+y(5)*y(2)+y(4)*y(3));
dydt2=[dy1;dy2;dy3;dy4;dy5;dy6;dy7];
end

%(rlvlh(1,:).*(y(4)^2-y(7)^2-y(6)^2+y(5)^2)+rlvlh(2,:).*(2*y(5)*y(6)-2*y(4)*y(7))+rlvlh(3,:).*2*y(4)*y(6))
%(rlvlh(1,:).*(2*y(6)*y(5)-2*y(4)*y(7))+rlvlh(2,:).*(y(6)^2-y(7)^2+y(4)^2-y(5)^2)+rlvlh(3,:).*(2*y(4)*y(5)-2*y(6)*y(7)))
%(rlvlh(1,:).*2*y(4)*y(6)+rlvlh(2,:).*(2*y(6)*y(7)-2*y(4)*y(5))+rlvlh(3,:).*(-q(7)^2-q(6)^2+q(5)^2+q(4)^2))