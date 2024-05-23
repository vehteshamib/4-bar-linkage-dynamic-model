clc
syms L1 L2 L3 L4 m1 m2 m3 T F P g real
syms q1 q2 q3 t real
syms u du real

q = [q1;q2;q3];

a =[ L1*cos(q1), L2*cos(q2), -L3*cos(q3)
     L1*sin(q1), L2*sin(q2), -L3*sin(q3)];
YY=[-1 1 0];

W =simplify(inv([YY;a])*[1;0;0])
dq= W*u;
dq1 = dq(1);  dq2 = dq(2);  dq3 = dq(3);
dW = zeros(3,1);
for i=1:3
   dW = dW +diff(W,q(i))*dq(i);
end

M =[      (L1^2*(m1 + 3*m2))/3, (L1*L2*m2*cos(q1 - q2))/2,           0
     (L1*L2*m2*cos(q1 - q2))/2,               (L2^2*m2)/3,           0
                             0,                         0, (L3^2*m3)/3]; 
B =[T - F*L1*cos(q1 - q2) - (L1*g*m1*sin(q1))/2 - L1*g*m2*sin(q1) - (L1*L2*dq2^2*m2*sin(q1 - q2))/2
                                                   -(L2*m2*(- L1*sin(q1 - q2)*dq1^2 + g*sin(q2)))/2
                                                                        (L3*(2*P - g*m3*sin(q3)))/2];
M_ML = simplify(W.'*M*W)
B_ML = simplify(W.'*(B-M*dW*u))







