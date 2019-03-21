function Q = qheart( t )
%Calculates an approximal function for LV inflow, just for testing!
Ts1=0.26;
Ts2=0.39;

for n=1:length(t)
    if mod(t(n),0.8)<Ts1
    Q(n)=1/2*(1-cos(mod(t(n),0.8)/Ts1*pi));
    end
    if mod(t(n),0.8)>=Ts1&&mod(t(n),0.8)<Ts2
    Q(n)=1/2*(1+cos((mod(t(n),0.8)-Ts1)/(Ts2-Ts1)*pi));   
    end
    if mod(t(n),0.8)>=Ts2
    Q(n)=0.001;
    end
end
Q=Q*500;
end

