function ea = eAtrium( t )
%Calculates an approximal function for LV inflow, just for testing!
Ts1=0.04;
Ts2=0.09;

    if mod(t,0.8)<Ts1
    ea=1/2*(1-cos(mod(t,0.8)/Ts1*pi));
    end
    if mod(t,0.8)>=Ts1&&mod(t,0.8)<Ts2
    ea=1/2*(1+cos((mod(t,0.8)-Ts1)/(Ts2-Ts1)*pi));   
    end
    if mod(t,0.8)>=Ts2
    ea=0;
    end


end

