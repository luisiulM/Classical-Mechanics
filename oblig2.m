clear all

v_r = 5;
A = 2;
b = 0.1;

phi = 0:0.0001:50;

V_r = @(phi) (v_r-b*phi)+A*cos(phi);

plot(phi, V_r(phi))
ylabel('radielt hastighet')
xlabel('fase')

%%
clear all

k = 8; % kg/s^2
m = 0.5; % kg

% Jeg velger beta slik at under resonansen, vil omega_d være reelt.
beta = 2; % kg/s
%

omega = k/m; % 1/s^2
gamma = beta/(2*m); % 1/s

% verdier til omega_d
omega_d = [1 6 20];
omega_d1 = sqrt( omega.^2 - 2*gamma.^2 );
%


t = 0:0.0001:10;

F = 1/m;
R = sqrt( (omega.^2-omega_d.^2).^2 + (2*gamma*omega_d).^2 );
R1 = sqrt( (omega.^2-omega_d1.^2).^2 + (2*gamma*omega_d1).^2 );
phi = atan( (2*gamma*omega_d)/(omega.^2-omega_d.^2) );

y = @(t2,R2,omega_d2) F/R2*cos(omega_d2*t - phi);

plot(t, y(t,R(1),omega_d(1)))
hold on
plot(t, y(t,R(2),omega_d(2)))
plot(t, y(t,R(3),omega_d(3)))
plot(t, y(t,R1,omega_d1))
xlabel('Tid')
ylabel('partikulær løsning y_p(t)')
legend('omega_d = 1', 'omega_d = 6', 'omega_d = 20', 'omega_d = sqrt( omega.^2 - 2*gamma.^2 )')

%%
clear all

m = 100;             % kg
R = 1.1*6.371*10.^6; % meter

% Her ganger jeg v med 60^2 for å gjøre om til timer
v = 6000*60^2;         % m/h

phi = pi/6;          % radianer
M = 5.972*10.^24;    % kg jordmasse

% Her ganger jeg G med 60^4 for å gjøre om til timer
G = 6.674*10.^-11*60.^4;   % m^3/kg*h^2
%

L = m*v*R*sin(phi);
a = G*M*m;
E = (1/2)*m*(v*cos(phi)).^2 + L.^2/(2*m*R.^2) - a/R;

delta_t = 0.00001;
t = 0:delta_t:5; % timer

theta = linspace(0,2*pi,1000);

% Start betingelser
r(1) = R;
O(1) = 0;
%

for i = 1:length(t)-1
    
    dr(i) = sqrt( 2/m*(E - L.^2/(2*m*r(i).^2) + a/r(i)) );
    dO(i) = L/(m*r(i).^2);
    
    r(i+1) = r(i) + dr(i)*delta_t;
    O(i+1) = O(i) + dO(i)*delta_t;
    
    if r(i) < 9.665*10.^6
        
       tid = t(i);
        
    end
    
end

% max radius
r_max = max(r);
%

% tid payloaden tar for å komme til r max
tid_h = tid;          % timer
tid_min = tid_h*60;   % minutter
tid_sec = tid_min*60; % sekunder
%


polarplot(O,r)
hold on

R2 = 6.371*10.^6;
r2 = R2*cos(theta).^2 + R2*sin(theta).^2;
polarplot(theta, r2)

R3 = 9.665*10.^6;
r3 = R3*cos(theta).^2 + R3*sin(theta).^2;
polarplot(theta, r3)
legend('payload', 'Jorda', 'payload på d)')
