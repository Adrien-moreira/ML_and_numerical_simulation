function [x,y,t] = RK4_2D(x0,y0,tmin,tmax,pas,F,G)
    t = tmin:pas:tmax;
    x = zeros(1,length(t));
    x(1) = x0;
    y = zeros(1,length(t));
    y(1) = y0;
    for k=1:length(t)-1
        k1f = F(t(k),x(k),y(k));
        k1g = G(t(k),x(k),y(k));
        k2f = F(t(k)+pas/2, x(k)+pas*k1f/2, y(k)+pas*k1g/2);
        k2g = G(t(k)+pas/2, x(k)+pas*k1f/2, y(k)+pas*k1g/2);
        k3f = F(t(k)+pas/2, x(k)+pas*k2f/2, y(k)+pas*k2g/2);
        k3g = G(t(k)+pas/2, x(k)+pas*k2f/2, y(k)+pas*k2g/2);
        k4f = F(t(k)+pas, x(k)+pas*k3f, y(k)+pas*k3g);
        k4g = G(t(k)+pas, x(k)+pas*k3f, y(k)+pas*k3g);
        
        % Mise à jour des coordonnees
        x(k+1) = x(k) + pas/6 * (k1f + 2*k2f + 2*k3f + k4f);
        y(k+1) = y(k) + pas/6 * (k1g + 2*k2g + 2*k3g + k4g);
    end
end
