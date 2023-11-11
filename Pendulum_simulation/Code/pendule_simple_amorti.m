clear variables;
close all;

eta = 1.8E-2;
r = 0.02;
m = 0.1;
l = 1;
g_pes = 9.81;

tmin = 0;
tmax = 150;

gamma = 6*pi*r*eta/m;
omega = sqrt(g_pes/l);

% Conditions initiales

theta_0 = 1;     % Position initiale
thetap_0 = 0;   % Vitesse initiale

f = @(t,theta,z)(z);
g = @(t,theta,z)(-gamma * z - omega^2 * sin(theta));

h = 0.07;   % Pas temporel

% 1. Méthode d'Euler
[thetaEuler,zEuler,tEuler] = Euler_2D(theta_0,thetap_0,tmin,tmax,h,f,g);

yEuler = -l*cos(thetaEuler);
xEuler = l*sin(thetaEuler);

% 2. Méthode RK2
betaRK = 1;
[thetaRK2,zRK2,tRK2] = RK2_2D(theta_0,thetap_0,tmin,tmax,h,betaRK,f,g);

yRK2 = -l*cos(thetaRK2);
xRK2 = l*sin(thetaRK2);

% 3. Méthode RK4
[thetaRK4,zRK4,tRK4] = RK4_2D(theta_0,thetap_0,tmin,tmax,h,f,g);

yRK4 = -l*cos(thetaRK4);
xRK4 = l*sin(thetaRK4);

% Animation
figure(1)
subplot(131);
H = plot(0, 0, 'Marker', 'o');
axis([-1.5 1.5 -1.5 1.5]);
axis equal;
title("Méthode d'Euler");
set(gca, 'nextplot', 'replacechildren');
grid on;

subplot(132);
H2 = plot(0, 0, 'Marker', 'o');
axis([-1.5 1.5 -1.5 1.5]);
axis equal;
title("Méthode RK2");
set(gca, 'nextplot', 'replacechildren');
grid on;

subplot(133);
H3 = plot(0, 0, 'Marker', 'o');
axis([-1.5 1.5 -1.5 1.5]);
axis equal;
title("Méthode RK4");
set(gca, 'nextplot', 'replacechildren');
grid on;

for i=1:min(length(tEuler),length(tRK2))
    
    % Pour terminer l'animation avant la fin de la boucle for
    
    if ~ishghandle(1)
        break;
    end
    
    set(H, 'Xdata', [0, xEuler(i)], 'Ydata',  [0, yEuler(i)]);
    set(H2, 'Xdata', [0, xRK2(i)], 'Ydata',  [0, yRK2(i)]);
    set(H3, 'Xdata', [0, xRK4(i)], 'Ydata',  [0, yRK4(i)]);
    drawnow;
end


% Visualisation de l'espace des phases

figure(2)
subplot(131);
plot(thetaEuler, zEuler);
xlabel("Theta");
ylabel("Theta prime");
title("Méthode d'Euler");
grid on;

subplot(132);
plot(thetaRK2, zRK2);
xlabel("Theta");
ylabel("Theta prime");
title("Méthode RK2");
grid on;

subplot(133);
plot(thetaRK4, zRK4);
xlabel("Theta");
ylabel("Theta prime");
title("Méthode RK4");
grid on;

% Calcul de l'énergie cinétique, potentielle et totale en fonction du temps

Ec_Euler = 0.5 * m * l^2 * zEuler.^2;
Ep_Euler = m * g_pes * l * (1 - cos(thetaEuler));
Etot_Euler = Ec_Euler + Ep_Euler;

Ec_RK2 = 0.5 * m * l^2 * zRK2.^2;
Ep_RK2 = m * g_pes * l * (1 - cos(thetaRK2));
Etot_RK2 = Ec_RK2 + Ep_RK2;

Ec_RK4 = 0.5 * m * l^2 * zRK4.^2;
Ep_RK4 = m * g_pes * l * (1 - cos(thetaRK4));
Etot_RK4 = Ec_RK4 + Ep_RK4;

% Énergie cinétique, potentielle et totale méthode d'Euler

figure(3)

subplot(1, 3, 1);
plot(tEuler, Ec_Euler);
xlabel("Temps (s)");
ylabel("E_c (J)");
title("Méthode d'Euler");

subplot(1, 3, 2);
plot(tEuler, Ep_Euler);
xlabel("Temps (s)");
ylabel("E_p (J)");
title("Méthode d'Euler");

subplot(1, 3, 3);
plot(tEuler, Etot_Euler);
xlabel("Temps (s)");
ylabel("E_{tot} (J)");
title("Méthode d'Euler");

% Énergie cinétique, potentielle et totale méthode RK2

figure(4)

subplot(1, 3, 1);
plot(tRK2, Ec_RK2);
xlabel("Temps (s)");
ylabel("E_c (J)");
title("Méthode RK2");

subplot(1, 3, 2);
plot(tRK2, Ep_RK2);
xlabel("Temps (s)");
ylabel("E_p (J)");
title("Méthode RK2");

subplot(1, 3, 3);
plot(tRK2, Etot_RK2);
xlabel("Temps (s)");
ylabel("E_{tot} (J)");
title("Méthode RK2");

% Énergie cinétique, potentielle et totale méthode RK4

figure(5)

subplot(1, 3, 1);
plot(tRK4, Ec_RK4);
xlabel("Temps (s)");
ylabel("E_c (J)");
title("Méthode RK4");

subplot(1, 3, 2);
plot(tRK4, Ep_RK4);
xlabel("Temps (s)");
ylabel("E_p (J)");
title("Méthode RK4");

subplot(1, 3, 3);
plot(tRK4, Etot_RK4);
xlabel("Temps (s)");
ylabel("E_{tot} (J)");
title("Méthode RK4");
