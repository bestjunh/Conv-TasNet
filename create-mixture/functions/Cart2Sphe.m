function [r, azi, ele] = Cart2Sphe(v)
x = v(1);
y = v(2);
z = v(3);

r = sqrt(x^2 + y^2 + z^2);
ele = acosd(z/r);
azi = atan2d(y, x);

end