clear all; clc; close all;

load('Rsdata4.mat');
load('calibration.mat');

M_PI = 3.141592653589793238;

acc.x_raw = rt_sensor.signals.values(:,1);
acc.y_raw = rt_sensor.signals.values(:,2);
acc.z_raw = rt_sensor.signals.values(:,3);
acc.time = rt_sensor.time;

gyro.p_raw = rt_sensor.signals.values(:,4);
gyro.q_raw = rt_sensor.signals.values(:,5);
gyro.r_raw = rt_sensor.signals.values(:,6);

%% Low Pass Filter

%Varyanslar üzerinden gürültü hesabı
%Değerlerin en stabil olduğu aralıklardan varyans hesabı yapıldı

acc.noise.x = var(acc.x_raw(3400:3850)); 
acc.noise.y = var(acc.y_raw(1600:3100));   
acc.noise.z = var(acc.z_raw(3400:4100));

gyro.noise.p = var(gyro.p_raw(1560:3120)); 
gyro.noise.q = var(gyro.q_raw(2100:3100)); 
gyro.noise.r = var(gyro.r_raw(1480:3120));

freq = 1/(acc.time(2)-acc.time(1));

%....................................................
acc.x = lowpass(acc.x_raw, acc.noise.x, freq);  
acc.y = lowpass(acc.y_raw, acc.noise.y, freq);  
acc.z = lowpass(acc.z_raw, acc.noise.z, freq);  

gyro.p = lowpass(gyro.p_raw, gyro.noise.p, freq);
gyro.q = lowpass(gyro.q_raw, gyro.noise.q, freq);
gyro.r = lowpass(gyro.r_raw, gyro.noise.r, freq);

%% Calibration
acc.measurement = [acc.x,acc.y,acc.z,ones(5600,1)];
acc.calibrated = acc.measurement*X_ls;

acc.x = acc.calibrated(:,1);
acc.y = acc.calibrated(:,2);
acc.z = acc.calibrated(:,3);

%% Raw Roll-Pitch values

for i=1:length(acc.y_raw)
roll_raw(i) = atan(acc.y_raw(i)/acc.z_raw(i));
end
roll_raw = roll_raw';

pitch_raw = (asin(acc.x_raw/9.80665)*180)/M_PI;

%% Kalman Filtresi
Va = 0;
g = 9.80665;
N = 5600;
Tout = 0.005;

phi = atan(acc.y(1)/acc.z(1));
theta = asin(acc.x(1)/g);

phi_measurement = rt_estim.signals.values(:,6);
theta_measurement = rt_estim.signals.values(:,5);

vari.phi = 10^-4;       %var(phi_measurement(950:1800))
vari.theta =  10^-4;     %var(theta_measurement(2800:3650))

R = [acc.noise.x 0 0; 0 acc.noise.y 0; 0 0 acc.noise.z];

x.previous = [phi,theta]';
P.previous = [vari.phi 0; 0 vari.theta];

Qu = [gyro.noise.p 0 0 0 ;
       0 gyro.noise.q 0 0;
       0 0  gyro.noise.r 0
       0 0 0 0];
   
Q_angle = [vari.phi 0; 0 vari.theta];
G = [1 sin(phi)*tan(theta) cos(phi)*tan(theta) 0;0 cos(phi) -sin(phi) 0];
Q = G*Qu*G' + Q_angle;

I = eye(2);
x.values(:,1) = x.previous;

% İlerletme (Prediction)

for i=2:N
    
    f = [gyro.p(i)+ gyro.q(i)*sin(phi)*tan(theta)+ gyro.r(i)*cos(phi)*tan(theta); 
                              gyro.q(i)*cos(phi)-gyro.r(i)*sin(phi)];
                                                        
    f_jacobian = [gyro.q(i)*cos(phi)*tan(theta)-gyro.r(i)*sin(phi)*tan(theta) , (gyro.q(i)*sin(phi)+gyro.r(i)*cos(phi))/(cos(theta))^2 ;
                                -gyro.q(i)*sin(phi)-gyro.r(i)*cos(phi),                                   0];

    x.predicted = x.previous + (Tout)*f;
    A = f_jacobian;
    P.predicted = P.previous + Tout*(A*P.previous + P.previous*A' + Q);
    
              
    h = [            gyro.q(i)*Va*sin(theta) + g*sin(theta);
         gyro.r(i)*Va*cos(theta)-gyro.p(i)*Va*sin(theta)-g*cos(theta)*sin(phi)
                     -gyro.q(i)*Va*cos(theta)-g*cos(theta)*cos(phi)];
    
    %Update state(Correction)
    
    x.current = x.predicted;
    P.current = P.predicted;
    
    if mod(i,1) == 0

    C = [0                                        (gyro.q(i)*Va*cos(theta)+g*cos(theta));
         (-g*cos(phi)*cos(theta))      (-gyro.r(i)*Va*sin(theta)-gyro.p(i)*Va*cos(theta)+g*sin(phi)*sin(theta));
         (g*sin(phi)*cos(theta))                 (gyro.q(i)*Va + g*cos(phi))*sin(theta)];
     
     
    L = (P.predicted*C')/(C*P.predicted*C'+R);
    %L = zeros(2,3);                              %İf only predicted state
    P.current = (I-L*C)*P.predicted;
     
    Y = [acc.x(i) acc.y(i) acc.z(i)]';
    
    x.current = x.predicted + L*(Y-h);  %Update
    
    end
    
    P.previous = P.current;
    
    x.values(:,i) = x.current;
    
    x.previous = x.current;
    P.previous = P.current;

    phi = x.current(1);
    theta = x.current(2);
    
end

x.val = x.values';

%% Plotting outputs

figure(1);
hold on;
a1 = plot(rad2deg(x.val(:,1))); 
M1 = "Kalman Çıktısı";
a2 = plot(rad2deg(roll_raw));
M2 = "Ham veri";
hold on;
legend([a1,a2], [M1, M2]);
title('Roll Açısı')
grid on;
xlabel('Zaman') 
ylabel('Derece °')

figure(2);
hold on;
a1 = plot(rad2deg(x.val(:,2))); 
M1 = "Kalman Çıktısı";
a2 = plot(pitch_raw);
M2 = "Ham veri";
hold on;
legend([a1,a2], [M1, M2]);
title('Pitch Açısı')
grid on;
xlabel('Zaman') ;
ylabel('Derece °');
