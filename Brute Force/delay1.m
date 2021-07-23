function [delay1] = delay1(QQ,SS,G,C)
x = (QQ.*C)./(SS.*G);
delay1=(.5*C*((1-(G/C)).^2))./(1-(min(1,x).*(G/C)));
end

