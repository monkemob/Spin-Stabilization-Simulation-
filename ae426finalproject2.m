clc
clear
close all
%same I as 1 but different M=[1;0;2]
time=[0:0.01:0.05];
time2=[0:.01:5];
[T,X5]=ode45(@diffeq5,time,[0;5;0;1;0;0;0]);
[T6,X6]=ode45(@diffeq6,time2,[X5(end,1);X5(end,2);5;X5(end,4);X5(end,5);X5(end,6);X5(end,7)]);
[T7,X7]=ode45(@diffeq7,time,[0;5;0;1;0;0;0]);
[T8,X8]=ode45(@diffeq8,time2,[X7(end,1);X7(end,2);5;X7(end,4);X7(end,5);X7(end,6);X7(end,7)]);

figure() 
plot(T6,X6(:,1))
title('Solution w');
xlabel('Time (s)');
ylabel('w');
hold on 
plot(T6,X6(:,2))
hold on
plot(T6,X6(:,3));
legend('wx','wy','wz')
title('body to lvlh')
hold off

figure()
plot(T6,X6(:,4))
hold on
plot(T6,X6(:,5))
hold on
plot(T6,X6(:,6))
hold on 
plot(T6,X6(:,7))
title('q WRT time')
xlabel('time (s)')
ylabel('q')
legend('q0','q1','q2','q3')
title('body to lvlh')

figure() 
plot(T8,X8(:,1))
title('Solution w');
xlabel('Time (s)');
ylabel('w');
hold on 
plot(T8,X8(:,2))
hold on
plot(T8,X8(:,3));
legend('wx','wy','wz')
title('body to eci')
hold off

figure()
plot(T8,X8(:,4))
hold on
plot(T8,X8(:,5))
hold on
plot(T8,X8(:,6))
hold on 
plot(T8,X8(:,7))
title('q WRT time')
xlabel('time (s)')
ylabel('q')
legend('q0','q1','q2','q3')
title('body to eci')

yourMatrix=[T7, X7(:,4), X7(:,5), X7(:,6), X7(:,7), zeros(6,1)+6600, zeros(6,1), zeros(6,1);T8, X8(:,4), X8(:,5), X8(:,6), X8(:,7),zeros(501,1)+6600, zeros(501,1), zeros(501,1) ];
writematrix(yourMatrix, 'output.csv', 'Delimiter', 'comma')

function dydt5=diffeq5(t,y) 
dy1=1-2*y(2)*y(3);
dy2=1.5*y(1)*y(3);
dy3=2-.25*y(1)*y(2);
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
dydt5=[dy1;dy2;dy3;dy4;dy5;dy6;dy7];
end

function dydt6=diffeq6(t,y)
dy1=-2*y(2)*y(3);
dy2=1.5*y(1)*y(3);
dy3=-.25*y(1)*y(2);
dy4=.5*(-y(5)*(y(1)-2*.00113*(y(5)*y(6)-y(4)-y(7)))-y(6)*(y(2)-.00113*(y(4)^2-y(5)^2+y(6)^2-y(7)^2))-y(7)*(y(3)-2*.00113*(y(6)*y(7)+y(4)*y(5))));
dy5=.5*(y(4)*(y(1)-2*.00113*(y(5)*y(6)-y(4)-y(7)))-y(7)*(y(2)-.00113*(y(4)^2-y(5)^2+y(6)^2-y(7)^2))+y(6)*(y(3)-2*.00113*(y(6)*y(7)+y(4)*y(5))));
dy6=.5*(y(7)*(y(1)-2*.00113*(y(5)*y(6)-y(4)-y(7)))+y(4)*(y(2)-.00113*(y(4)^2-y(5)^2+y(6)^2-y(7)^2))-y(5)*(y(3)-2*.00113*(y(6)*y(7)+y(4)*y(5))));
dy7=.5*(-y(6)*(y(1)-2*.00113*(y(5)*y(6)-y(4)-y(7)))+y(5)*(y(2)-.00113*(y(4)^2-y(5)^2+y(6)^2-y(7)^2))+y(4)*(y(3)-2*.00113*(y(6)*y(7)+y(4)*y(5))));
dydt6=[dy1;dy2;dy3;dy4;dy5;dy6;dy7];
end

function dydt7=diffeq7(t,y) 
%angular velocity
dy1=1-2*y(2)*y(3);
dy2=1.5*y(1)*y(3);
dy3=2-.25*y(1)*y(2);
%qs. this format is same for all diffeq functions
dy4=.5*(-y(5)*y(1)-y(6)*y(2)-y(7)*y(3));
dy5=.5*(y(4)*y(1)-y(7)*y(2)+y(6)*y(3));
dy6=.5*(y(7)*y(1)+y(4)*y(2)-y(5)*y(3));
dy7=.5*(-y(6)*y(1)+y(5)*y(2)+y(4)*y(3));
dydt7=[dy1;dy2;dy3;dy4;dy5;dy6;dy7];
end
%torque free
function dydt8=diffeq8(t,y) 
dy1=-2*y(2)*y(3);
dy2=1.5*y(1)*y(3);
dy3=-.25*y(1)*y(2);
dy4=.5*(-y(5)*y(1)-y(6)*y(2)-y(7)*y(3));
dy5=.5*(y(4)*y(1)-y(7)*y(2)+y(6)*y(3));
dy6=.5*(y(7)*y(1)+y(4)*y(2)-y(5)*y(3));
dy7=.5*(-y(6)*y(1)+y(5)*y(2)+y(4)*y(3));
dydt8=[dy1;dy2;dy3;dy4;dy5;dy6;dy7];
end
