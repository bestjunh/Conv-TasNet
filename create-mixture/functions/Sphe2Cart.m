function v = Sphe2Cart(dist,azi,ele)

x = dist*sind(ele)*cosd(azi);
y = dist*sind(ele)*sind(azi);
z = dist*cosd(ele);

v = [x y z].';

end