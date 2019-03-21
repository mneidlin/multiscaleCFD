function ev = eVentricle( t )
%Calculates an approximal function for LV inflow, just for testing!
Ts1=0.26;
Ts2=0.39;

    if mod(t,0.8)<Ts1
    ev=1/2*(1-cos(mod(t,0.8)/Ts1*pi));
    end
    if mod(t,0.8)>=Ts1&&mod(t,0.8)<Ts2
    ev=1/2*(1+cos((mod(t,0.8)-Ts1)/(Ts2-Ts1)*pi));   
    end
    if mod(t,0.8)>=Ts2
    ev=0;
    end


end

