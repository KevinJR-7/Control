clc
clear 
close all

%% Solución por Matlab, punto 5

syms s t I1 I2 I3
syms E R L C1 C2 positive

% Ecuaciones de Impedancias en el Circuito
E1 = (L*s + 1/(s*C1))*I1 - (1/(C1*s))*I2 - (L*s)*I3 - E/s == 0;
E2 = -(1/(C1*s))*I1 + (1/(C1*s) + 2*R)*I2 - (R)*I3 == 0;
E3 = -L*s*I1 - R*I2 + (L*s + 1/(s*C2) + R)*I3 == 0;

% Resolver el sistema de ecuaciones
sol = solve([E1, E2, E3], [I1, I2, I3]);

% Calcular la Inversa de Laplace de las corriente
I1_t = ilaplace(sol.I1);
I2_t = ilaplace(sol.I2);
I3_t = ilaplace(sol.I3);
%%

% Valores Definidos para las constantes
E_val = 5;
L_val = 1;
C1_val = 10e-6;
C2_val = 10e-6;
R_val = 1000;

%%
% Sustituir valores en la expresión de I(t)
I1_t_expr = subs(I1_t, [E, L, C1, C2, R], [E_val, L_val, C1_val, C2_val, R_val]);
I2_t_expr = subs(I2_t, [E, L, C1, C2, R], [E_val, L_val, C1_val, C2_val, R_val]);
I3_t_expr = subs(I3_t, [E, L, C1, C2, R], [E_val, L_val, C1_val, C2_val, R_val]);

%%
% Crear un vector de tiempo
t_values = linspace(0, 0.1, 1000); % Definir un rango de tiempo

% Crear un Vector de Ceros
I1_values = zeros(size(t_values));
I2_values = zeros(size(t_values));
I3_values = zeros(size(t_values));

% Evaluar las expresiones en el Vector de Tiempos
for i = 1:length(t_values)    
    I1_values(i) = subs(I1_t_expr, t, t_values(i));
    I2_values(i) = subs(I2_t_expr, t, t_values(i));
    I3_values(i) = subs(I3_t_expr, t, t_values(i));
end
%% Hallar de Forma numérica el voltaje en los Capacitores y el Inductor

Vc1 = (1/C1_val) * cumtrapz(t_values, I1_values - I2_values);
Vc2 = (1/C2_val) * cumtrapz(t_values, I3_values);
Vl = L_val*(diff(I1_values - I3_values ) ./ diff(t_values));

%%

% Crear la figura y los subplots
figure;

% Primer subplot
subplot(3, 2, 1); 
plot(t_values, I1_values, 'b', 'LineWidth', 2);
xlabel('Tiempo');
ylabel('I_{1}(t)');
title('I_{1}(t)');
grid on;

% Segundo subplot
subplot(3, 2, 2); 
plot(t_values, Vc1, 'r', 'LineWidth', 2);
xlabel('Tiempo');
ylabel('V_{C1}(t)');
title('V_{C1}(t)');
grid on;

% Tercer subplot
subplot(3, 2, 3); 
plot(t_values, I2_values, 'g', 'LineWidth', 2);
xlabel('Tiempo');
ylabel('I_{3}(t)');
title('I_{2}(t)');
grid on;

% Cuarto subplot
subplot(3, 2, 4); 
plot(t_values, Vc2, 'm', 'LineWidth', 2);
xlabel('Tiempo');
ylabel('V_{C2}(t)');
title('V_{C2}(t)');
grid on;

% Quinto subplot
subplot(3, 2, 5); 
plot(t_values, I3_values, 'c', 'LineWidth', 2);
xlabel('Tiempo');
ylabel('I_{3}(t)');
title('I_{3}(t)');
grid on;

% Sexto subplot
subplot(3, 2, 6); 
plot(t_values(1:end-1), Vl, 'k', 'LineWidth', 2);
xlabel('Tiempo');
ylabel('V_{L}(t)');
title('V_{L}(t)');
grid on;
