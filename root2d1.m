function G = root2d1(x)
global L1 L2 L3 L4 q1_0 dq1_0
q1 = q1_0;
dq1 = dq1_0;
aaa = x(1);
G(1)= L1*sin(q1) - L4 + L2*sin(aaa) - L3*sin(x(2));
G(2)= L3*cos(x(2)) - L2*cos(x(1)) - L1*cos(q1);
G(3)= L1*dq1*cos(q1) + L2*x(3)*cos(x(1)) - L3*x(4)*cos(x(2));
G(4)= L1*dq1*sin(q1) + L2*x(3)*sin(x(1)) - L3*x(4)*sin(x(2));