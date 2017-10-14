%% This is a script that simulates the motion of a quadcopter, using a 
%% PID controller and a sliding mode controller

%% Parameters and states initializing
clc; clear;

% parameters

m = 4.5;		% mass in kilogram
g = 9.8;	% acceleration of gravity
Ixx = 0.05;	% 
Iyy = 0.05;
Izz = 0.1;
kdx = 0.8;
kdy = 0.8;
kdz = 1;
Kd = diag([kdx,kdy,kdz]);
I = diag([Ixx,Iyy,Izz]);
s = zeros(4,1);

% states
x = 0;
y = 0;
z = 0;
phi = 0.5;
theta = -0.4;
psi = 1;
u = 0;
v = 0;
w = 0;
p = 0;
q = 0;
r = 0;

X1 = [x;y;z];
X2 = [u;v;w];
X3 = [phi;theta;psi];
X4 = [p;q;r];

start_time = 0;
end_time = 20;
dt = 0.01;
timePeriod = start_time:dt:end_time;

fid1 = fopen('states.txt','w');
fid2 = fopen('input.txt','w');

% Desired states

X3d = [0;0;0];
e1 = 0;
inte1 = 0;
currentez = X1(3);
intez = 0;

%% Quadcopter states updating by equation

for t = timePeriod
    if t == 2.33
        
    end
    phi = X3(1);
    theta = X3(2);
    psi = X3(3);
	R_omega = [1 0 -sin(theta); 0 cos(phi) cos(theta)*sin(theta); 0 -sin(phi) cos(phi) * cos(theta)];	% rotation matrix of angular velocity
	R = [cos(theta)*cos(psi)-cos(theta)*sin(phi)*sin(psi) -cos(psi)*sin(phi) - cos(phi)*cos(theta)*sin(psi) sin(theta)*sin(psi);...
         cos(theta)*cos(psi)*sin(phi) + cos(phi)*sin(psi) cos(phi)*cos(theta)*cos(psi) - sin(phi)*sin(psi) -cos(psi)*sin(theta);...
         sin(phi)*cos(theta) cos(phi)*sin(theta) cos(theta)];		% rotation matrix of coordinates
	interMediate1 = [((Iyy-Izz)/Ixx)*X4(2)*X4(3);((Izz-Ixx)/Iyy)*X4(1)*X4(3);((Ixx-Iyy)/Izz)*X4(1)*X4(2)];

% 	Determine control input
        
	laste1 = e1;
	e1 = X3d - X3;
	inte1 = inte1 + dt*e1;
	de1 = (e1-laste1)/dt;
    
	%ControlInput = PIDcontrol(t, e1, de1, inte1, end_time, X3, m, g, I);
    %ControlInput(2:4) = sat(ControlInput(2:4),0.5);
	ControlInput = SMControl(X1,X2,X3,X4, interMediate1, R_omega, m, g, kdz, I);

%	States update
    mdisturb = 0;
    %idisturb = [0.5*sin(t), 0.4*cos(t), 0.4*sin(t+0.3)];
    
	dX1 = X2;
	dX2 = [0;0;-g] + (R/(m + mdisturb)) * [0;0;ControlInput(1)] - Kd*X2*(1/(m + mdisturb)) ;
	dX3 = inv(R_omega) * X4 ;
	dX4 = inv(I) * ControlInput(2:4) - interMediate1 + 4*sin(3*t)*[1;1;1];

	X1 = X1 + dt*dX1;
	X2 = X2 + dt*dX2;
	X3 = X3 + dt*dX3;
	X4 = X4 + dt*dX4;

% Data processing
	%fprintf(fid1,'%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n',t,X1,X2,X3,X4)
	fprintf(fid1,'%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n',t,X1,X3);
    fprintf(fid2,'%8.4f %8.4f %8.4f %8.4f %8.4f\n',t,ControlInput);
end

DATA = load('states.txt');
DATA2 = load('input.txt');

figure(1);
plot(DATA(:,1),DATA(:,5),'r'); title('phi angle');xlabel('time/s');ylabel('angle/rad');ylim([-1 1]);hold on; plot(DATA(:,1),zeros(numel(DATA(:,1)),1),'b--');
%print phi -dpng;
figure(2);
plot(DATA(:,1),DATA(:,6),'r'); title('theta angle');xlabel('time/s');ylabel('angle/rad');ylim([-1 1]);hold on; plot(DATA(:,1),zeros(numel(DATA(:,1)),1),'b--');
%print theta -dpng;
figure(3);
plot(DATA(:,1),DATA(:,7),'r'); title('psi angle');xlabel('time/s');ylabel('angle/rad');ylim([-1 1]);hold on; plot(DATA(:,1),zeros(numel(DATA(:,1)),1),'b--');
%print psi -dpng;


figure(4);
plot(DATA2(:,1),DATA2(:,2),'b'); title('Z input');xlabel('time/s');ylabel('Drag force/N');ylim([40 50]);
%print Thrust -dpng;
figure(5);
plot(DATA2(:,1),DATA2(:,3),'b'); title('X moment');xlabel('time/s');ylabel('Inertia/kg*m^2');
%print Xmoment -dpng;
figure(6);
plot(DATA2(:,1),DATA2(:,4),'b'); title('Y moment');xlabel('time/s');ylabel('Inertia/kg*m^2');
%print Ymoment -dpng;
figure(7);
plot(DATA2(:,1),DATA2(:,5),'b'); title('Z moment');xlabel('time/s');ylabel('Inertia/kg*m^2');
%print Zmoment -dpng;



function output = PIDcontrol(t,e1,de1,inte1,intermediate1,X3,m,g,I)	% Generate a PID control scheme
	kp = diag([3 3 3]);
	ki = diag([0.01 0.01 0.01]);
	kd = diag([2 2 2]);
    output = zeros(4,1);
    output(1) = m*g/(cos(X3(1))*cos(X3(2)));

	if t >= 0
		output(2:4) = kp*e1 + kd*de1 + ki*inte1;
	else
		output(2:4) = kp*e1 + kd*de1;
    end
end

function output = SMControl(X1,X2,X3,X4, intermediate1, R_omega, m, g, kdz, I)	% Generate a sliding mode scheme
	k1 = 1; k2 = 1; k3 = 1; k4 = 1;
    beta1 = 2; beta2 = 2; beta3 = 1.5;
	output  = zeros(4,1);
    s = zeros(4,1);

	s(1) = X2(3) + k1*X1(3);
	s(2) = X4(1) + k2*X3(1);
	s(3) = X4(2) + k3*X3(2);
	s(4) = X4(3) + k4*X3(3);

	%output(1) = (m/cos(X3(2)) * (g + (kdz-k1)*X2(3))) - beta1*sign(s(1));
    output(1) = m*g/(cos(X3(1))*cos(X3(2)));
	output(2:4) = I * (intermediate1 - R_omega * diag([k2,k3,k4]) * X3) - diag([beta1,beta2,beta3])*sat(s(2:4), 10);
end	

function res = sat(x, A)
    res = x;
    for i = 1:numel(x)
        res(i) = 2/(1+exp(-A*x(i))) - 1;
    end
end