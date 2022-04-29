clear all ; close all;  clc;

% Given Properties
a    = 0;       b    = 0.1;
h    = 0.001;   yINI = 0;
t    = a:h:b;   y(1) = yINI;
N    = (b-a)/h;

tau  = 0.01;    T    = 1/tau;
f    = 100;     Vm   = 1;
w    = 2*pi*f;

% Solving Equation Using Ode 45
[tmat,ymat] = ode45(@myRC, [a b], yINI); % Fourth/Fifth RK
figure();
plot(tmat,ymat,'*k');
hold on;

% Euler Explicit Method
yE(1)      = yINI;
for      i = 1:N
    t(i+1) = t(i) + h;
   yE(i+1) = yE(i) + myRC(t(i),yE(i))*h;
end
plot(t,yE,'-b');
hold on;

% Euler Modified Method
yEu    = 0;
Slope2 = 0;
yEM(1) = yINI;
for      i = 1:N
    t(i+1) = t(i) + h;
    Slope1 = myRC(t(i),yEM(i));
       yEu = yEM(i)+Slope1*h;
    Slope2 = myRC(t(i+1), yEu);
  yEM(i+1) = yEM(i) + (Slope1+Slope2)*h/2;
end
plot(t,yEM,'-r')
hold off

% Plot Propertiews
title("NM ODE Euler");
grid on;
xlabel('x'); ylabel('y');
xlim([0 0.1])
ylim([-0.02, 0.02])
legend('ODE45' , 'EU' , 'EM');

% Given Problem
function dvdx = myRC(t,v)
 tau  = 1;
 T    = 1/tau;
 f    = 10;
 Vm   = 1;
    dvdx = -T*v + 1*T*Vm*cos(2*pi*f*t);
end
