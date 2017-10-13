%% This is a script that simulates the motion of a quadcopter, using a 
%% PID controller and a sliding mode controller

%% Parameters and states initializing

global I, kdx, kdy, kdz, Kd, Ixx, Iyy, Izz, interMediate1

% parameters

m = 4.5;		% mass in kilogram
g = 9.8;	% acceleration of gravity
Ixx = 3;	% 
Iyy = 3;
Izz = 3;
kdx = 1.7;
kdy = 0.8;
kdz = 0.65;
Kd = diag([kdx,kdy,kdz]);
I = diag([Ixx,Iyy,Izz]);
beta1 = 1;
beta2 = 1;
beta3 = 1;
beta4 = 1;
s = zero(4,1);

% states
x = 0;
y = 0;
z = 0;
phi = 0;
theta = 0;
psi = 0;
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
end_time = 10;
dt = 0.01;
timePeriod = start_time:dt:end_time;

fid1 = fopen('states.txt');

% Desired states

X3d = [0;0;0];

%% Quadcopter states updating by equation

for t = timePeriod
	R_omega = [1;2;3];	% rotation matrix of angular velocity
	R = [1;2;3];		% rotation matrix of coordinates
	interMediate1 = [((Iyy-Izz)/Ixx)*X4(2)*X4(3);((Izz-Ixx)/Iyy)*X4(2)*X4(4);((Ixx-Iyy)/Izz)*X4(3)*X4(4)];

% 	Determine control input
	laste1 = e1;
	e1 = X3d - X3;
	inte1 = inte1 + dt*e1;
	de1 = (e1-laste1)/dt;
	Controlinput = PIDcontrol(t, e1, de1, inte1);
	Controlinput = SMControl(X1, X2, X3, X4);

%	States update

	dX1 = X2;
	dX2 = [0;0;-g] + R * [0;0;U1] - Kd*X2 ;
	dX3 = inv(R_omega) * X4 ;
	dX4 = inv(I) * ControlInput - interMediate1;

	X1 = X1 + dt*dX1;
	X2 = X2 + dt*dX2;
	X3 = X3 + dt*dX3;
	X4 = X4 + dt*dX4;

% Data processing
	%fprintf(fid1,'%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n',t,X1,X2,X3,X4)
	fprintf(fid1,'%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n',t,X1,X2);
end

function output = PIDcontrol(t,e1,de1,inte1)	% Generate a PID control scheme
	kp = 3;
	ki = 5.5;
	kd = 4;

	if t >= end_time/2
		output = kp*e1 + kd*de1 + ki*inte1;
	else
		output = kp*e1 + kd*de1;
	end

end

function output = SMControl(X1,X2,X3,X4)	% Generate a sliding mode scheme
	k1 = 1; k2 = 1; k3 = 1; k4 = 1;
	output  = zeros(4,1);

	s(1) = X2(3) + k1*X1(3);
	s(2) = X4(1) + k2*X3(1);
	s(3) = X4(2) + k3*X3(2);
	s(4) = X4(3) + k4*X3(3);

	output(1) = (m/cos(X3(2)) * (g + (kdz-k1)*X2(3)) - beta1*sign(s(1))
	output(2:4) = I * (intermediate1 - R_omega * diag([k2,k3,k4]) * X3) - diag([beta1,beta2,beta3])*sign(s(2:4));

	
end	