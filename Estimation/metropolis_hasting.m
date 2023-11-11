function metropolis_hasting(T, a, b, sigma_q)
    
    x = zeros(1, T); 
    x_cand = zeros(1, T); 
    x(1) = 0.5; 
    
    % Calcul du rapport d'acceptation 
    for i = 2:T
        x_cand(i) = x(i-1) + sigma_q * randn();
    
        if x_cand(i) < 0 || x_cand(i) > 1
            alpha = 0; 
        else
            ratio_proposition = 1; 
            alpha = min(1, (betapdf(x_cand(i), a, b) / betapdf(x(i-1), a, b)) * ratio_proposition);
        end
    
        if rand() < alpha
            x(i) = x_cand(i);       
        else
            x(i) = x(i-1); 
        end 
    end 
    
    % Estimation de la densité de probabilité empirique
    [count, center] = hist(x, 100); 
    dx = center(2) - center(1); 
    bar(center, count / (T * dx)); 
    hold on;
    
    % Tracé de la loi cible théorique
    fplot(@(x) betapdf(x, a, b), [0 1], 'r');
    legend('Estimée', 'Théorique');
    title('Densité de probabilité');
    hold off;

    % Tracé de la série temporelle des x
    figure;
    plot(1:T, x);
    title('Série Temporelle');
    
    % Auto-corrélation
    figure;
    [param1, param2] = xcorr(x-mean(x), 'coeff');
    plot(param2, param1);
    title('Fonction d''auto-corrélation');
end
