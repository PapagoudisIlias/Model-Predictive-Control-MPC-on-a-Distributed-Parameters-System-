clear all
% close all
clc

%Make sure to add the path to Casadi correctly on your computer
addpath('C:\Users\ippap\Documents\CASADI\casadi-windows-matlabR2014b-v3.5.5')
import casadi.*

N=10;
Nc=5;
j=3;
NN=3;
mm=2^j+1;

dx=1/2^j;
dT=0.01;
mpc_dT=10*dT;
simtime=3;

%Matrices phiINT and DphiINT are essential for the space discretization
phiINT=zeros(2^j+1,2^j);
for i=2:9
    phiINT(i,i-1)=1;
end

DphiINT=zeros(2^j+1,2^j);
DphiINT(1,:)=[-3.72453,1.02455,-0.258108,0.045612,-0.00533224,0.000405926,-0.0000200259,6.39115*10^(-7)];
DphiINT(2,:)=[-1.64445, -6.04501, 1.58732, -0.354709, 0.0566837, -0.00616414,0.000446618, -0.0000213173];
DphiINT(3,:)=[6.76679, -0.057131, -6.39971, 1.644, -0.360873, 0.0571303,-0.00618546, 0.000447283];
DphiINT(4,:)=[-1.70158, 6.41208, -0.000447284, -6.40588, 1.64445, -0.360894, 0.057131, -0.00618547];
DphiINT(5,:)=[0.360894, -1.64445, 6.4059, 0, -6.4059, 1.64445, -0.360894, 0.057131];
DphiINT(6,:)=[-0.057131, 0.360894, -1.64445, 6.40588, 0.000447284, -6.41208, 1.70158, -0.721788];
DphiINT(7,:)=[0.00618546, -0.0571303, 0.360873, -1.644, 6.39971, 0.057131,-6.76679, 3.2889];
DphiINT(8,:)=[-0.000446618, 0.00616414, -0.0566837, 0.354709, -1.58732, 6.04501, 1.64445, -12.8118];
DphiINT(9,:)=[0.0000200259, -0.000405926, 0.00533224, -0.045612, 0.258108,-1.02455, 3.72453, 10.2447];




sta = SX.sym('sta',2^j+1,1); 
states = [sta]; n_states=length(states);

v = SX.sym('v'); 
controls = [v]; n_controls = length(controls);


rhs=[0];

for k=2:(2^j+1)
    sumphi=0;
    for m=0:2^j
        sumphi=sumphi+sta(m+1)*phiINT(m+1,k-1);
    end
    Dsumphi=0;
    for m=0:2^j
        Dsumphi=Dsumphi+v*sta(m+1)*DphiINT(m+1,k-1);
    end
    rhs=[rhs;-sumphi-Dsumphi];
end

f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U',n_controls,Nc); % Decision variables (controls)
P = SX.sym('P',n_states + n_states);

X = SX.sym('X',n_states,(N+1));

% compute solution symbolically

X(:,1) = P(1:mm); % initial state
for k = 1:N
    st = X(:,k); 
    if k<= Nc
        con = U(:,k);
    else
        con=U(:,Nc); % constsnt input after Nc timesteps
    end
    
    for cc = 1:(mpc_dT/dT)
        k1=f(st,con);
        k2=f(st+k1*(dT/2),con);
        k3=f(st+k2*(dT/2),con);
        k4=f(st+k3*dT,con);
        st=st+ dT*(k1+2*k2+2*k3+k4)/6;
    end
    st_next = st;
    X(:,k+1) = st_next;
end

ff=Function('ff',{U,P},{X});

obj = 0; % Objective function
g = [];  % constraints vector

qd=zeros(1,mm); qd(1)=1; qd((mm-1)/4+1)=1; qd((mm-1)/2+1)=10; qd(3*(mm-1)/4+1)=50;  qd(mm)=50;
Q=diag([qd]);
R=[10];


%Objective function
for k=1:N
    st = X(:,k);  
    if k<=Nc-1
        con = U(:,k);
        con_next = U(:,k+1);
        obj = obj+(st-P(n_states+1:2*n_states))'*Q*(st-P(n_states+1:2*n_states)) + (con-con_next)'*R*(con-con_next); % calculate obj
    else
        obj = obj+(st-P(n_states+1:2*n_states))'*Q*(st-P(n_states+1:2*n_states)) ; % After Nc timesteps con=constant so con-con_next=0
    end
end

% compute constraints
for k = 1:N+1   
    for p=1:mm
        g=[g;X(p,k)];
    end
end

OPT_variables = U;
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;

% constraints of control
args.lbx = 0.9;
args.ubx = 1.1;


args.lbg=-10;
args.ubg=10;

t0=0;
x0 = ones(n_states,1);


xs=zeros(n_states,1);


%Reference temperatures

xs(1)=1;xs((mm-1)/4+1)=0.778801; xs((mm-1)/2+1)=0.606531; xs(3*(mm-1)/4+1)=0.472367;xs(mm)=0.367879;



xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = ones(Nc,1);

sim_tim = 10; % Maximum simulation time


% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];

xtime=zeros(1,n_states);
time=zeros(1,n_states);

xlength=[];
for i=0:n_states-1
    xlength=[xlength i*dx];
end

ramplength=xlength;
xtemp=ones(1,n_states);


sigma=0*0.008;

maxreps=(1/mpc_dT)*simtime+1;
main_loop = tic;

while (mpciter < sim_tim / dT && mpciter<maxreps)
    args.p   = [x0;xs]; % set the values of the parameters vector
    args.x0 = u0;
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);    
    u=full(sol.x);
    
    %In case artificial noise in control action is needed
%     ra=normrnd(0,sigma);
%     u(1)=u(1)+ra;
%     if u(1)<0.9 | u(1)>1.1
%         u(1)=u(1)-ra;
%     end

 %In case the first control must be u=1
%     if mpciter==0
%         u(1)=1;
%     end
    
    ff_value = ff(u',args.p);
    xx1(:,1:n_states,mpciter+1)= full(ff_value)';
    
    
    t(mpciter+1) = t0;
    % In case Δu constraints are needed (Δu = u(t+1)-u(t))
%     if mpciter==0
%         uprev=1;
%     else
%         uprev=u_cl(end);
%     end
%     dumaxmin=0.005;
%     if (u(1)-uprev>dumaxmin)  
%         u(1)=uprev+dumaxmin;
%     elseif (u(1)-uprev<-dumaxmin)
%         u(1)=uprev-dumaxmin;
%     end
%     
    u_cl= [u_cl ; u(1)];
    
    [t0, x0, u0] = shift_DPS_bigdt_RK4(mpc_dT,dT, t0, x0, u,f,sigma);
    
    %In case artificial noise in measurements is needed we add white noise
%     x0((mm-1)/4+1)=x0((mm-1)/4+1)+normrnd(0,sigma);
%     x0((mm-1)/2+1)=x0((mm-1)/2+1)+normrnd(0,sigma);
%     x0(3*(mm-1)/4+1)=x0(3*(mm-1)/4+1)+normrnd(0,sigma);
%     x0(mm)=x0(mm)+normrnd(0,sigma);



    
    xx(:,mpciter+2) = x0(:,1);  
    
    xtime=xtime+dT;
    time=[time xtime];
    ramplength=[ramplength xlength];
    xtemp=[xtemp x0(:,1)'];
    
    disp(['mpciter = ',num2str(mpciter)])
    mpciter = mpciter + 1;
end
main_loop_time = toc(main_loop)
avg_loop_time=main_loop_time/maxreps
disp(' ')
    output=norm([x0(1);x0((mm-1)/4+1);x0((mm-1)/2+1);x0(3*(mm-1)/4+1);x0(mm)]-[xs(1);xs((mm-1)/4+1);xs((mm-1)/2+1);xs(3*(mm-1)/4+1);xs(mm)]);
    refer=norm([xs(1);xs((mm-1)/4+1);xs((mm-1)/2+1);xs(3*(mm-1)/4+1);xs(mm)]);
    disp(['The error is ',num2str(output/refer*100),' %'])
    [ [x0(1);x0((mm-1)/4+1);x0((mm-1)/2+1);x0(3*(mm-1)/4+1);x0(mm)],[xs(1);xs((mm-1)/4+1);xs((mm-1)/2+1);xs(3*(mm-1)/4+1);xs(mm)]]
figure(1)
t=0:mpc_dT:maxreps*mpc_dT-mpc_dT;
stairs(t,u_cl')
title(['Control actions'])
xlim([0 t(end)])
xlabel('time')
ylabel('velocity')
grid on

figure(2)
[lenx leny]=size(time);
for n=1:leny/n_states-1
    clf
    plot(ramplength(1:n_states),xtemp((n-1)*n_states+1:(n*n_states)),'-*')
    text(0.1,0.45,['time = ' num2str((n-1)*mpc_dT)])
    hold on
    scatter([0.25 0.5 0.75 1],[xs((mm-1)/4+1) xs((mm-1)/2+1) xs(3*(mm-1)/4+1) xs(mm)],'o','r')
    hold off
    grid on
    title(['Real time output'])
    xlabel('ramplength')
    ylabel('Temperature')
    ylim([0 1])

    movieVector(n) = getframe;
    pause(0.1)
end


xtemp(end-mm+1:end);

figure(3)
xlen=0:1/8:1;
plot(xlen,xtemp(end-mm+1:end),'-*')
hold on
scatter([0.25 0.5 0.75 1],[xs((mm-1)/4+1) xs((mm-1)/2+1) xs(3*(mm-1)/4+1) xs(mm)],'o','r')
title(['Final output'])
xlabel('ramplength')
ylabel('temperature')
xlim([0,1.2])
ylim([0.2,1])
grid on
