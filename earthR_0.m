function R_0=earthR_0(A,B,x)
%发射系下地心矢径,A为发射方位角，B为纬度,x为模式，1为输入地理纬度，2为输入地心纬度，3为球形模型
    ae=6578140;
    be=6356755;
    r0=6.37111e6;
    if x==1
    a=ae*be;
    phi=atand(be^2/ae^2 * tand(B));
    b=sqrt(be^2 * cosd(phi)^2 + ae^2 * sind(phi)^2);
    r_0=a/b;%地心矢长度
    miu=B-phi;
    R_0=r_0*[-sind(miu)*cosd(A);
             cosd(miu);
             sind(miu)*sind(A)];
    elseif x==2
    a=ae*be;
    phi=atand(ae^2/be^2 * tand(B));
    b=sqrt(be^2 * cosd(B)^2 + ae^2 * sind(B)^2);
    r_0=a/b;%地心矢长度
    miu=phi-B;
    R_0=r_0*[-sind(miu)*cosd(A);
             cosd(miu);
             sind(miu)*sind(A)];
    elseif x==3
    R_0=r0*[0,1,0];
    end
end