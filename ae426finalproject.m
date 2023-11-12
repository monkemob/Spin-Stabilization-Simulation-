clc
clear
close all

time=[0:0.01:0.1];
time2=[0:.05:15];
[T,X]=ode45(@diffeq,time,[0;0;5;1;0;0;0]);
[T2,X2]=ode45(@diffeq2,time2,[X(end,1);X(end,2);5;X(end,4);X(end,5);X(end,6);X(end,7)]);
[T3,X3]=ode45(@diffeq3,time,[0;0;5;1;0;0;0]);
[T4,X4]=ode45(@diffeq4,time2,[X3(end,1);X3(end,2);5;X3(end,4);X3(end,5);X3(end,6);X3(end,7)]);
yourMatrix=[T4, X4(:,4), X4(:,5), X4(:,6), X4(:,7),zeros(301,1)+6600, zeros(301,1), zeros(301,1) ];

writematrix(yourMatrix, 'output2.csv', 'Delimiter', 'comma')

%precession and spin rate
Hg=[100*X4(:,1)'; 200*X4(:,2)'; 400*X4(:,3)';];
theta=[1:301];
for k=1:301
    thetavec(k)=acos(dot(Hg(:,k),[0,0,1])/norm(Hg(:,k)));
end
%since ode45 is numerical, it is not exact so there is some variation in
%theta
theta=mean(thetavec);
%same as above for omega
omegavec=X4(:,1)./cos(9.*T4);
omega=mean(omegavec);
%spin rate
psidot=5-omega/tan(theta);
%precession rate
phidot=400/(100-400)*psidot/cos(theta);

figure() 
plot(T2,X2(:,1))
title('Solution w');
xlabel('Time (s)');
ylabel('w');
hold on 
plot(T2,X2(:,2))
hold on
plot(T2,X2(:,3));
legend('wx','wy','wz')
title('body to lvlh')
hold off

figure()
plot(T2,X2(:,4))
hold on
plot(T2,X2(:,5))
hold on
plot(T2,X2(:,6))
hold on 
plot(T2,X2(:,7))
title('q WRT time')
xlabel('time (s)')
ylabel('q')
legend('q0','q1','q2','q3')
title('body to lvlh')

figure() 
plot(T4,X4(:,1))
title('Solution w');
xlabel('Time (s)');
ylabel('w');
hold on 
plot(T4,X4(:,2))
hold on
plot(T4,X4(:,3));
legend('wx','wy','wz')
title('body to eci')
hold off

figure()
plot(T4,X4(:,4))
hold on
plot(T4,X4(:,5))
hold on
plot(T4,X4(:,6))
hold on 
plot(T4,X4(:,7))
title('q WRT time')
xlabel('time (s)')
ylabel('q')
legend('q0','q1','q2','q3')
title('body to eci')

%distrubance
function dydt=diffeq(t,y) 
dy1=1-2*y(2)*y(3);
dy2=2+1.5*y(1)*y(3);
dy3=-.25*y(1)*y(2);
dy4=.5*(-y(5)*(y(1)-2*.00113*(y(5)*y(6)-y(4)-y(7)))-y(6)*(y(2)-.00113*(y(4)^2-y(5)^2+y(6)^2-y(7)^2))-y(7)*(y(3)-2*.00113*(y(6)*y(7)+y(4)*y(5))));
dy5=.5*(y(4)*(y(1)-2*.00113*(y(5)*y(6)-y(4)-y(7)))-y(7)*(y(2)-.00113*(y(4)^2-y(5)^2+y(6)^2-y(7)^2))+y(6)*(y(3)-2*.00113*(y(6)*y(7)+y(4)*y(5))));
dy6=.5*(y(7)*(y(1)-2*.00113*(y(5)*y(6)-y(4)-y(7)))+y(4)*(y(2)-.00113*(y(4)^2-y(5)^2+y(6)^2-y(7)^2))-y(5)*(y(3)-2*.00113*(y(6)*y(7)+y(4)*y(5))));
dy7=.5*(-y(6)*(y(1)-2*.00113*(y(5)*y(6)-y(4)-y(7)))+y(5)*(y(2)-.00113*(y(4)^2-y(5)^2+y(6)^2-y(7)^2))+y(4)*(y(3)-2*.00113*(y(6)*y(7)+y(4)*y(5))));
% y(8)=2*.00113*(y(5)*y(6)-y(4)-y(7));
% y(9)=.00113*(y(4)^2-y(5)^2+y(6)^2-y(7)^2);
% y(10)=2*.00113*(y(6)*y(7)+y(4)*y(5));
% y(11)=y(1)-y(8);
% y(12)=y(2)-y(9);
% y(13)=y(3)-y(10);
dydt=[dy1;dy2;dy3;dy4;dy5;dy6;dy7];
end
%torque free
function dydt2=diffeq2(t,y)
dy1=-2*y(2)*y(3);
dy2=1.5*y(1)*y(3);
dy3=-.25*y(1)*y(2);
dy4=.5*(-y(5)*(y(1)-2*.00113*(y(5)*y(6)-y(4)-y(7)))-y(6)*(y(2)-.00113*(y(4)^2-y(5)^2+y(6)^2-y(7)^2))-y(7)*(y(3)-2*.00113*(y(6)*y(7)+y(4)*y(5))));
dy5=.5*(y(4)*(y(1)-2*.00113*(y(5)*y(6)-y(4)-y(7)))-y(7)*(y(2)-.00113*(y(4)^2-y(5)^2+y(6)^2-y(7)^2))+y(6)*(y(3)-2*.00113*(y(6)*y(7)+y(4)*y(5))));
dy6=.5*(y(7)*(y(1)-2*.00113*(y(5)*y(6)-y(4)-y(7)))+y(4)*(y(2)-.00113*(y(4)^2-y(5)^2+y(6)^2-y(7)^2))-y(5)*(y(3)-2*.00113*(y(6)*y(7)+y(4)*y(5))));
dy7=.5*(-y(6)*(y(1)-2*.00113*(y(5)*y(6)-y(4)-y(7)))+y(5)*(y(2)-.00113*(y(4)^2-y(5)^2+y(6)^2-y(7)^2))+y(4)*(y(3)-2*.00113*(y(6)*y(7)+y(4)*y(5))));
dydt2=[dy1;dy2;dy3;dy4;dy5;dy6;dy7];
end

%body to eci. since wr=0 wb/r=wb 
function dydt3=diffeq3(t,y) 
dy1=1-2*y(2)*y(3);
dy2=2+1.5*y(1)*y(3);
dy3=-.25*y(1)*y(2);
dy4=.5*(-y(5)*y(1)-y(6)*y(2)-y(7)*y(3));
dy5=.5*(y(4)*y(1)-y(7)*y(2)+y(6)*y(3));
dy6=.5*(y(7)*y(1)+y(4)*y(2)-y(5)*y(3));
dy7=.5*(-y(6)*y(1)+y(5)*y(2)+y(4)*y(3));
dydt3=[dy1;dy2;dy3;dy4;dy5;dy6;dy7];
end
%torque free
function dydt4=diffeq4(t,y) 
dy1=-2*y(2)*y(3);
dy2=+1.5*y(1)*y(3);
dy3=-.25*y(1)*y(2);
dy4=.5*(-y(5)*y(1)-y(6)*y(2)-y(7)*y(3));
dy5=.5*(y(4)*y(1)-y(7)*y(2)+y(6)*y(3));
dy6=.5*(y(7)*y(1)+y(4)*y(2)-y(5)*y(3));
dy7=.5*(-y(6)*y(1)+y(5)*y(2)+y(4)*y(3));
dydt4=[dy1;dy2;dy3;dy4;dy5;dy6;dy7];
end
