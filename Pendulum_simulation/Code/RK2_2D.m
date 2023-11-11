function [x,y,t] = RK2_2D(x0,y0,tmin,tmax,pas,beta,F,G)
    t = tmin:pas:tmax;
    x = zeros(1,length(t));
    x(1) = x0;
    y = zeros(1,length(t));
    y(1) = y0;
    for k=1:length(t)-1
        k1f = F(t(k),x(k),y(k));
        k1g = G(t(k),x(k),y(k));
        k2f = F(t(k)+pas/(2*beta),x(k)+pas*k1f/(2*beta),y(k)+pas*k1g/(2*beta));
        k2g = G(t(k)+pas/(2*beta),x(k)+pas*k1f/(2*beta),y(k)+pas*k1g/(2*beta));
        
        % Mise à jour des coordonnees
        x(k+1) = x(k)+pas*((1-beta)*k1f+beta*k2f);
        y(k+1) = y(k)+pas*((1-beta)*k1g+beta*k2g);
    end
end