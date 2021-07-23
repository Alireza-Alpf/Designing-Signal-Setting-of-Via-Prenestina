function [delay2] = delay2(QQ,SS,G,C,T)
x = (QQ.*C)./(SS.*G);
c=G.*SS./C;
delay2 = (900*T*((x-1)+sqrt(((x-1).^2)+((4*x)./(c*T)))));
end 