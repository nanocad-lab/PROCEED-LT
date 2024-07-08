function leakage=Leakage_fun(b,Vdd)

leakage=b(1).*Vdd.^6+b(2).*Vdd.^5+b(3).*Vdd.^4+b(4).*Vdd.^3+b(5).*Vdd.^2+b(6).*Vdd+b(7);