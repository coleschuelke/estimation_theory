function plotcar_cole(x, h)
% h is now the plot handle

%calculate the car corners
Rot = [cos(x(4)), -sin(x(4)); ...
       sin(x(4)), cos(x(4))];
%the corners of the car: (front/back) (driver/passenger)
fd = x(1:2) + Rot*[x(5); x(6)/2];
fp = x(1:2) + Rot*[x(5); -x(6)/2];
bd = x(1:2) + Rot*[0; x(6)/2];
bp = x(1:2) + Rot*[0; -x(6)/2];

h.XData = [fp(1), fd(1), bd(1), bp(1), fp(1)];
h.YData = [fp(2), fd(2), bd(2), bp(2), fp(2)];

return;
