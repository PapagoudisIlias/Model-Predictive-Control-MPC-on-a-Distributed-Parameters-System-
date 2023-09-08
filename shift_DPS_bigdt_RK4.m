function [t0,x0,u0]=shift_DPS_bigdt_RK4(mpc_dT,dT,t0,x0,u,f,sigma)
st=x0;
con=u(1);


for cc=1:(mpc_dT/dT)
   k1=f(st,con);
   k2=f(st+k1*(dT/2),con);
   k3=f(st+k2*(dT/2),con);
   k4=f(st+k3*dT,con);
   st=st+ dT*(k1+2*k2+2*k3+k4)/6;
   % disturbances in initial velocity
   con=u(1)+normrnd(0,sigma);
end

x0=full(st);

t0=t0+dT;
u0=[u(:,2:end),u(end)]';
end
