function A=goaround(a,axis)
%a为转动角度，axis为转动轴,1,2,3代表x,y,z轴
if axis == 1
    A=[1, 0, 0;
       0, cosd(a), sind(a);
       0, -sind(a), cosd(a)];
elseif axis == 2
    A=[cosd(a), 0, -sind(a);
       0, 1, 0;
       sind(a), 0, cosd(a)];
elseif axis == 3
    A=[cosd(a), sind(a), 0;
       -sind(a), cosd(a), 0;
       0, 0, 1];
end
end