function avgDelay=TDi(C,g,Q,S,T)
         
         x = (Q.*C)./(S.*g);
         d1=(.5*C*((1-(g/C)).^2))./(1-(min(1,x).*(g/C)));
         c=g.*S./C;
         d2 = (900*T*((x-1)+sqrt(((x-1).^2)+((4*x)./(c*T)))));
         
         avgDelay=(d1+d2);
        
end