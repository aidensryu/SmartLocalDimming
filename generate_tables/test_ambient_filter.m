% Test temporal filter

freq = 10; % Sampling freqency [Hz]

cutoff_freq = 0.5; %Hz

N = 30*freq;
E_s = 10.^(rand(N,1)*2);
rep = round(rand(N,1)*20);

E = zeros(size(E_s));
t = linspace(0, N*1/freq, N );

prev_kk=1;
E(1)=E_s(1);
for kk=1:N
    if( kk>=(rep(prev_kk)+prev_kk) )
        E(kk) = E_s(kk);
        prev_kk = kk;
    else
        E(kk) = E(prev_kk);
    end
end

%%

c_f = 0.5/(freq/2);
f = [0 c_f c_f 1];
m = [1 1 0 0];
b = fir2(30,f,m);

c_f = 2/(freq/2);
f = [0 c_f c_f 1];
m = [1 1 0 0];
b_high = fir2(10,f,m);

%[b,a] = butter(30,0.5/(freq/2));

pad = length(b)-1;
pad_high = length(b_high)-1;
Ex = padarray( E, [pad 0], E(1), 'pre' );
Exf = zeros(size(E));

Exf(1) = E(1);
deltaT = 1000/60;
tau = 20;
for kk=2:length(E)
    T = E(kk);
    deltaE = (T-Exf(kk-1));
%    if( deltaE < 0 )
%        tau = 5;
%    else
%        tau = 7;
%    end
    Exf(kk) = T - deltaE*exp(-deltaT/tau);
end

Exf2 = zeros(length(E)/2,1);
Exf2(1) = E(1);
deltaT = 1000/30;
tau = 200;
for kk=2:1:(length(E)/2)
    T = E(kk*2);
    deltaE = (T-Exf2(kk-1));
    Exf2(kk) = T - deltaE*exp(-deltaT/tau);
end


Ef = filter(b,1,E, repmat(E(1),[length(b)-1,1]));

F = E.^0.5;
Ff = (max(filter(b,1,F, repmat(F(1),[length(b)-1,1])),1)).^2;

clf;

plot( t, E, '-b' )
hold on
plot( t, Ef, '-r' )
plot( t, Exf, '-g' )
plot( t(1:2:end), Exf2, '--m' )

hold off
set( gca, 'YScale', 'log' );