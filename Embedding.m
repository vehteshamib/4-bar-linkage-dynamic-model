function dz = Embedding (t,z)

global L1 L2 L3 L4 m1 m2 m3 T F P g
q1 = z(1);  q2 = z(2);  q3 = z(3);
u= z(4);

W =[                                    1
     -(L1*sin(q1 - q3))/(L2*sin(q2 - q3))
     -(L1*sin(q1 - q2))/(L3*sin(q2 - q3))];
 
 
M_ML = -(L1^2*(2*m1 + 5*m2 + 2*m3 - 3*m2*cos(2*q1 - 2*q2) - 2*m1*cos(2*q2 - 2*q3) + m2*cos(2*q1 - 2*q3) - 2*m3*cos(2*q1 - 2*q2) - 3*m2*cos(2*q2 - 2*q3)))/(6*(cos(2*q2 - 2*q3) - 1));
 
 
B_ML = -(6*F*L1*L2*L3*cos(q1 - q2)*sin(q2 - q3)^4 - 6*L2*L3*T*sin(q2 - q3)^4 - 2*L1^3*L2*m3*u^2*cos(q2 - q3)*sin(q1 - q2)^3 + 2*L1^3*L3*m2*u^2*cos(q2 - q3)*sin(q1 - q3)^3 + 6*L1*L2*L3*P*sin(q1 - q2)*sin(q2 - q3)^3 - 2*L1^3*L2*m2*u^2*cos(q2 - q3)*sin(q1 - q2)*sin(q1 - q3)^2 + 2*L1^3*L3*m3*u^2*cos(q2 - q3)*sin(q1 - q2)^2*sin(q1 - q3) + 3*L1*L2*L3*g*m1*sin(q2 - q3)^4*sin(q1) + 6*L1*L2*L3*g*m2*sin(q2 - q3)^4*sin(q1) + 3*L1^3*L3*m2*u^2*sin(q1 - q2)*sin(q1 - q3)^2*sin(q2 - q3)^2 - 3*L1^3*L2*m2*u^2*cos(q1 - q2)*cos(q1 - q3)*sin(q1 - q2)*sin(q2 - q3)^2 - 3*L1^3*L3*m2*u^2*cos(q1 - q2)*cos(q2 - q3)*sin(q1 - q3)^2*sin(q2 - q3) - 3*L1^2*L2*L3*m2*u^2*cos(q1 - q2)*cos(q1 - q3)*sin(q2 - q3)^3 + 2*L1^2*L2*L3*m3*u^2*cos(q1 - q2)*sin(q1 - q2)*sin(q2 - q3)^2 + 2*L1^2*L2*L3*m2*u^2*cos(q1 - q3)*sin(q1 - q3)*sin(q2 - q3)^2 + 3*L1^2*L2*L3*m2*u^2*sin(q1 - q2)*sin(q1 - q3)*sin(q2 - q3)^3 - 3*L1*L2*L3*g*m2*sin(q1 - q3)*sin(q2 - q3)^3*sin(q2) - 3*L1*L2*L3*g*m3*sin(q1 - q2)*sin(q2 - q3)^3*sin(q3) + 2*L1^3*L2*m2*u^2*cos(q1 - q3)*sin(q1 - q2)*sin(q1 - q3)*sin(q2 - q3) + 2*L1^3*L3*m3*u^2*cos(q1 - q2)*sin(q1 - q2)*sin(q1 - q3)*sin(q2 - q3) + 3*L1^3*L2*m2*u^2*cos(q1 - q2)*cos(q2 - q3)*sin(q1 - q2)*sin(q1 - q3)*sin(q2 - q3))/(6*L2*L3*sin(q2 - q3)^4);


X = M_ML\B_ML;

dz = [W*u ; X];