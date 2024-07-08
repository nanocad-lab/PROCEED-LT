function Vth = BTI(type, Vdd, Vth0, temp)
temp = temp + 275;
t = 1e8;
if type ==0
    gama = 5.2;
    alpha = 0.158;
    A = 3.12e-2;
    k = 50;
    Tin = 1.4;
else
    gama = 3;
    alpha = 0.173;
    A = 2.02e-2;
    k = 50;
    Tin = 1.4;
end
Eox = (Vdd - Vth0)/Tin;
Vth = A*exp(-k/temp)*(t^alpha)*(Eox^gama);
end