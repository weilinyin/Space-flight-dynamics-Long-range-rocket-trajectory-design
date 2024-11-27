function A=goaround(a,axis)
%a为转动角度，axis为转动轴,1,2,3代表x,y,z轴
if axis == 1
    A=[1, 0, 0;
       0, cos(a), sin(a);
       0, -sin(a), cos(a)];
elseif axis == 2
    A=[cos(a), 0, -sin(a);
       0, 1, 0;
       sin(a), 0, cos(a)];
elseif axis == 3
    A=[cos(a), sin(a), 0;
       -sin(a), cos(a), 0;
       0, 0, 1];
end
end