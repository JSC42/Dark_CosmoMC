function [r]=Interp_2D(y_axis_,x_axis_,Y,X,z,method)
% Interpolation for 2D array
% Written by Jason @ IHEP, Beijing
% Works well if size(z) > 10*10
% INPUT:
% -- x_axis,y_axis: x and y array
% -- z: z Array, elements given by z(y_index,x_index)
% -- [y,x]: where to interpolate
% -- Interp_METHOD:
% ---------- 1: Linear
% ---------- 2: Log Linear

FileExist=exist('Find_Index.m');
if ~FileExist
    ! cp ~/FUNCTIONS/Find_Index.m ./
end
if method ==1
    x = X;
    y = Y;
    y_axis_=y_axis_;
    x_axis_=x_axis_;
else
    x = log10(X);
    y = log10(Y);
    y_axis_=log10(y_axis_);
    x_axis_=log10(x_axis_);
end

xid1 = Find_Index(x_axis_,x);
yid1 = Find_Index(y_axis_,y);

xid2=xid1+1;
yid2=yid1+1;

x1=x_axis_(xid1);
x2=x_axis_(xid2);
y1=y_axis_(yid1);
y2=y_axis_(yid2);

f11=z(yid1,xid1);
f12=z(yid1,xid2);
f21=z(yid2,xid1);
f22=z(yid2,xid2);

f1=((f21-f11)/(y2-y1))*(y-y1)+f11;
f2=((f22-f12)/(y2-y1))*(y-y1)+f12;

f=((f2-f1)/(x2-x1))*(x-x1)+f1;
r=f;

end
