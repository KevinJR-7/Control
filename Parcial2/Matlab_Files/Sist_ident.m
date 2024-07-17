clc;
clear all;
close all;

%% Datos Depurados
data = readtable('Datos/Data_end.csv');

ur = data{:,"Channel2_V_"};
ur = ur(1:10:end);

out = data{:,"Channel1_V_"};
out = out(1:10:end);

time = data{:,"Time_s_"};
time = time(1:10:end);

Tm = time(2)-time(1);

%Otro modelo
% data = readtable('Datos/complete_df.csv');
% ur = data{:,"Channel1_V_"};
% ur = ur(1:20:end);
% out = data{:,"Channel2_V_"};
% out = out(1:20:end);
% time = data{:,"Time_s_"};
% time = time(1:20:end);


datos = iddata(out,ur,Tm); % Datos ident

%El modelo de LSM es un modelo muy sensible a los datos redundantes,
%por lo que es importante diezmar los datos es decir reducir la 
%frecuencia de muestreo

%% Gráfica de la Salida

figure,

hold on;

plot(time,ur,'LineWidth', 2, 'DisplayName', 'u(t)')
plot(time,out,'LineWidth', 2, 'DisplayName', 'y(t)')

hold off;

xlabel('Tiempo(s)');
ylabel('Voltaje(V)');
title('Output');
legend('show','Location', 'northeast');

grid on;

%% Matrices del sistema (Identificación LSM)

% Este método se utiliza para encontrar la función de transferencia del 
% sistema discreto, hallando los coeficientes de un modelo lineal que 
% representan la curva de mejor ajuste a los datos

%Esto se realiza modelando la función de transferencia en una matriz de 
%ecuaciones en Diferencias de la forma 

% y[k] = -a1*y[k-1] - a2*y[k-2] + b0*u[k-1] + b1*u[k-2]   

%Para un sistema de Orden 2 con máximo retraso 2 

% {y[k]}   = {-y[k-1] -y[k-2] u[k-1] u[k-2]} {a1}
% {y[k-1]} = {-y[k-2] -y[k-3] u[k-2] u[k-3]} {a2}
% {y[k-2]} = {-y[k-3] -y[k-4] u[k-3] u[k-4]} {b0}
% {y[k-3]} = {-y[k-4] -y[k-5] u[k-4] u[k-5]} {b1}

%Donde a son los coeficientes de la salida y b los de la entrada
%Y es el vector de entradas 
%PHI es la matriz conformada por los yk y uk (Observabilidad - Determ)
%theta es el vector de coeficientes

% Y se puede escribir de forma: Y = Phi*theta

% Despejando los coeficientes llegamos a la forma 
% theta = inv(Phi'*Phi)*(Phi'*Y); 

N = length(ur);
Y = out(3:N);
Phi = [-out(2:N-1) -out(1:N-2) ur(2:N-1) ur(1:N-2)];


%% Icógnitas

theta = inv(Phi'*Phi)*(Phi'*Y);

Gz_approx = tf(theta(3:4)',[1 theta(1:2)'],Tm);

num = cell2mat(Gz_approx.Numerator);
den = cell2mat(Gz_approx.Denominator);

[Adl,Bdl,Cdl,Ddl] = tf2ss(num, den);

[y1, t_d, x1] = lsim(Gz_approx, ur, time);

%% Con Ident
 
Options = tfestOptions;                        
Options.EnforceStability = true;               
                                            
tf_id = tfest(datos, 2, 1, Options, 'Ts', Tm);

%% Gráficos

[y2, t, x2] = lsim(tf_id, ur, time);

figure,
hold on;

plot(t, y1,'LineWidth', 2, 'DisplayName', 'LSM');
plot(t, y2,'LineWidth', 2, 'DisplayName', 'Ident');
plot(t, out,'LineWidth', 2, 'DisplayName', 'Real(t)');

hold off;

xlabel('Tiempo');
ylabel('Respuesta');
title('Respuesta al escalón');
legend('show')
grid on;

%% Métricas de Ajuste Identificación

Tabla_error(out,y1,'LSM')
Tabla_error(out,y2,'Ident')

%% Control por Realimentación de estados

ts = 16/2;
Mp = 0.28/2;

zita = sqrt(((log(Mp))^2)/(((log(Mp))^2)+ pi^2));
wn = 4.6/(ts*zita);


%Hallar los polos continuos deseados y discretizarlos
polos = roots([1 2*zita*wn wn^2]);
p_discretos = exp(polos*Tm);


%% Espacio de Estados Discreto

[Ad,Bd,Cd,Dd] = tf2ss(tf_id.Numerator, tf_id.Denominator);

ssd = ss(Ad,Bd,Cd,Dd,Tm);

Kd = place(Ad,Bd,p_discretos);

%Kd Ident
%Kd = place(Ad,Bd,p_discretos);

%% Metricas de Desempeño Reales

%Sistema Modelado con Retroalimentación
Act = Ad - Bd*Kd;
ssd = ss(Act,Bd,Cd,Dd,Tm);

%Datos de la tarjeta
DB = readtable('Datos/DB.csv');
Board = DB{:,1};

rd = step(ssd,time);

Tabla_error(Board,rd,'Adquisición')

%% Espacio de Estados Continuo

tsc = 16;
Mpc = 0.28;

%Polos reales
zitac = sqrt(((log(Mpc))^2)/(((log(Mpc))^2)+ pi^2));
wnc = 4.6/(tsc*zitac);

%Polos Deseados
zitac_d = sqrt(((log(Mpc/2))^2)/(((log(Mpc/2))^2)+ pi^2));
wnc_d = 4.6/(tsc/2*zitac_d);

pc = roots([1 2*zitac_d*wnc_d wnc_d^2]);

[A,B,C,D] = tf2ss(wnc^2, [1 2*zitac*wnc wnc^2]); 

Kdl = place(A,B,pc);

%% Función para calcular errores

function Tabla_error(pv_f,datos,Name)

    MAPE = (1/length(pv_f))*sum(abs((pv_f - datos)./pv_f)*100);
    [R2,RMSE] = rsquare(pv_f,datos);

    % Tabla con los resultados
    resultados = table({Name}, MAPE, R2, RMSE, 'VariableNames', {'Modelo', 'MAPE', 'R2', 'RMSE'});

    % Muestra la tabla
    disp('Resultados Estadísticos:');
    disp(resultados);
end
