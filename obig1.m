%                                  Oppgave 4

clear all

%m = 2; % kg
%k = 8; % N/m
%y0 = 0.5; % m
%g = 9.81; % m/s^2

%dt1 = 0.000001; % sec
%dt2 = 0.1; %sec

%t1 = 0:dt1:10;
%t2 = 0:dt2:10;

%y1(1) = 0.3; % m
%v1(1) = 0; % m/s
%a1(1) = -k/m * (y1(1) - y0) - g; % m/s^2

%for i = 1:length(t1)-1
    
   %y1(i+1) = y1(i) + v1(i) * dt1;
   %v1(i+1) = v1(i) + a1(i) * dt1;
   %a1(i+1) = -k/m * (y1(i) - y0) - g;
    
%end

%y2(1) = 0.3;
%v2(1) = 0;
%a2(1) = -k/m * (y2(1) - y0) - g;

%for i = 1:length(t2)-1
    
    %y2(i+1) = y2(i) + v2(i) * dt2;
    %v2(i+1) = v2(i) + a2(i) * dt2;
    %a2(i+1) = -k/m * (y2(i) - y0) - g;
    
%end
 
%subplot(1,3,1)
%plot(t1,y1)
%xlabel('Tid')
%ylabel('y(t)')
%title('y(t) med dt = 0.000001 sec')
%hold on

%subplot(1,3,2)
%plot(t2,y2)
%xlabel('Tid')
%ylabel('y(t)')
%title('y(t) med dt = 0.1 sec')

%A = 1.125;

%Y = 2*A*cos(sqrt((k/m))*t1);

%subplot(1,3,3)
%plot(t1,Y)
%xlabel('Tid')
%ylabel('Y(t)')
%title('Y(t)= 2*A*cos(sqrt((k/m)*t))')

%                                 Oppgave 1

%L1 = 2; % meter
%L2 = 6; % meter
%m = 1000; % kg
%g = 9.81; % m/s^2

%a1 = linspace(0,pi,50);

%S = @(a) L1*m*g*sin(a)/L2;

%plot(a1, S(a1))
%xlabel('vinkel')
%ylabel('Snordrag')


%                                 Oppgave 3

v0 = 9*sqrt(218);
y0 = 1;
g = 9.81;
o1 = linspace(0,pi/2,1000);

x = @(o) v0*cos(o) .* ( v0*sin(o) + sqrt( (v0*sin(o)).^2 + 2*g*y0 ))/g;

plot(o1,x(o1))
xlabel('vinkel')
ylabel('Rekkevidde')