% Mixotroph evolution model v2: temperature scaling for respiration

function dxdt = mixo_model(t,x,st,Tt,It,p,opt_time)

% define parameters

a  = p(1);
b1 = p(2);
b2 = p(3);
u1 = p(4);
u2 = p(5);
p1 = p(6);
p2 = p(7);
g1 = p(8);
g2 = p(9);
r1 = p(10);
r2 = p(11);
r3 = p(12);
lambda = p(13);
gam = p(14);
Qminxc = p(15);
Qminxn = p(16);
bc = p(17);
bn = p(18);
ki = p(19);
kb = p(20);
c = p(21);
T0 = p(22);

% interpolate light and temperature to model time

T = interp1(st, Tt, t);
I = interp1(st, It, t);

% variables

B    = x(1); 
Xa   = x(2);
Xp   = x(3);
Xg   = x(4);
That = x(5);
Qxc  = x(6);
Qxn  = x(7);
X    = x(8);

% update optimization strategy

global opt_R

nextR = find(isnan(opt_R(:,1)),1);
if isempty(nextR)
    nextR = length(opt_R);
end
if t >= opt_time(nextR)
    opt_R(nextR,:) = optimize_growth;
    R = opt_R(nextR,:);
else
    R = opt_R(nextR-1,:);
end

% adjust rates for temperature and investment strategy

ub = b1*b2^((T-T0)/10);
ux = Xa*u1*u2^((T-T0)/10);
px = Xp*p1*p2^((T-T0)/10);
gx = Xg*g1*g2^((T-T0)/10);
rx = r1*(T - That)^2 + r2*r3^((T-T0)/10);

% equations

dxdt = zeros(8,1);
dxdt(1) = ub*B - gx*X*B/(kb + B) - a*B;
dxdt(2) = (R(1) - Xa)*c; 
dxdt(3) = (R(2) - Xp)*c; 
dxdt(4) = (R(3) - Xg)*c;
dxdt(5) = (T - That)*lambda;
dxdt(6) = px*I/(ki + I) + gx*bc*B/(kb + B) - ux*min([1 - Qminxc/Qxc, 1 - Qminxn/Qxn],[],2)*Qxc;
dxdt(7) = gx*bn*B/(kb + B) + gam*rx*Qxn - ux*min([1 - Qminxc./Qxc, 1 - Qminxn./Qxn],[],2)*Qxn;
dxdt(8) = ux*min([1 - Qminxc/Qxc, 1 - Qminxn/Qxn],[],2)*X - rx*X - a*X;

% growth optimization function

function R_opt = optimize_growth

Ra = [];
Rp = [];
Rg = [];
syms Ra Rp Rg

ux = Ra*u1*u2^((T-T0)/10);
px = Rp*p1*p2^((T-T0)/10);
gx = Rg*g1*g2^((T-T0)/10);
rx =r1*(T - That)^2 + r2*r3^((T-T0)/10);
Qcstar = 1/ux*(px*I/(ki + I) + gx*bc*B/(kb + B) + ux*Qminxc);
Qnstar = 1/(ux - gam*rx)*(gx*bn*B/(kb + B) + ux*Qminxn);
mum_max = [ux*(1 - Qminxc/Qcstar); ux*(1 - Qminxn/Qnstar)];

fm = matlabFunction(-1*mum_max, 'Vars', {[Ra Rp Rg]});
x0 = [1/3 1/3 1/3];
options = optimoptions('fminimax','Display','off');
R_opt = fminimax(fm, x0, [], [], [1 1 1], 1, [0 0 0], [], [], options);

end


end