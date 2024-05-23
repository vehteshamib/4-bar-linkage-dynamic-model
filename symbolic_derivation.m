clc
syms L1 L2 L3 L4 m1 m2 m3 T F P g real
syms q1 q2 q3 t real
syms dq1 dq2 dq3 real

q = [q1;q2;q3];
dq = [dq1;dq2;dq3];

r_G1 = [L1/2*sin(q1);L1/2*-cos(q1)];
r_A  = [L1*sin(q1);L1*-cos(q1)];
r_G2 = [L1*sin(q1);L1*-cos(q1)] + [L2/2*sin(q2);L2/2*-cos(q2)];
r_B  = [L4;0] + [L3*sin(q3);L3*-cos(q3)];
r_G3 = [L4;0] + [L3/2*sin(q3);L3/2*-cos(q3)];

V_G1 = jacobian(r_G1,q)*dq;
V_A  = jacobian(r_A,q) *dq;
V_G2 = jacobian(r_G2,q)*dq;
V_B  = jacobian(r_B,q) *dq;
V_G3 = jacobian(r_G3,q)*dq;

% Generalized Forces
FF = F*[-cos(q2);-sin(q2)];
PP = P*[cos(q3);sin(q3)];
W = FF.'*V_A + PP.'*V_B + T*dq1;
Q = jacobian(W,dq).';

% Kinetic Energy
K = (1/2*m1*V_G1.'*V_G1 + 1/2*(1/12*m1*L1^2)*dq1^2) + (1/2*m2*V_G2.'*V_G2 + 1/2*(1/12*m2*L2^2)*dq2^2) + (1/2*m3*V_G3.'*V_G3 + 1/2*(1/12*m3*L3^2)*dq3^2);

% Potential Energy
U = m1*g*r_G1(2) + m2*g*r_G2(2) + m3*g*r_G3(2);

% Lgrangian
L=K-U;
E=K+U;

dL_ddq = jacobian(L,dq).';
M = jacobian(dL_ddq,dq);
B = Q - (jacobian(dL_ddq,q)*dq + diff(dL_ddq,t) - jacobian(L,q).');

M = simplify(M);
B = simplify(B);
E = simplify(E);

% Constraints
C = [L1*sin(q1);L1*-cos(q1)] + [L2*sin(q2);L2*-cos(q2)] - [L3*sin(q3);L3*-cos(q3)] - [L4;0];
dC = jacobian(C,q)*dq;

a = simplify(jacobian(dC,dq));
a_dot = zeros(2,3);
for i=1:3
    a_dot=a_dot+diff(a,q(i))*dq(i);
end
a_dot_dq = simplify(a_dot*dq);

% Embedding
W = [-inv(a(:,1:2))*a(:,3);1];
dW = jacobian(W,q)*dq;
M_Emb = W.'*M*W;
B_Emb = W.'*B - W.'*M*dW*q3;

W = simplify(W);
M_Emb = simplify(M_Emb);
B_Emb = simplify(B_Emb);

% Energy error
E = simplify(E);
W = simplify(Q.'*dq);
