
%Recursive Lagrange for standard Denavit-Hartenberg notation

function [G3] = gravity_L(q)


th1=q(1);
th2=q(2);
th3=q(3);
th4=q(4);
th5=q(5);
th6=q(6);
th7=q(7);

% standard DH Parameters of the robot
d = [0.2703 0 0.3644 0 0.3743 0 0.2295];
a = [0.069 0 0.069 0 0.01 0 0];
th = [th1 th2 th3 th4 th5 th6,th7]; 
al = [-pi/2 pi/2 -pi/2 pi/2 -pi/2 pi/2 0];


% Generalized transformation matrix for a link CS
% Tr = [cos(th) -sin(th)*cos(al)  sin(th)*sin(al) a*cos(th); ...
%       sin(th)  cos(th)*cos(al) -cos(th)*sin(al) a*sin(th); ...
%             0  sin(al)                  cos(al)         d; ...
%             0        0                       0          1]

%Forward kinematic transformation  matrices for individual link

A1 = [cos(th(1))  -sin(th(1))*cos(al(1))  sin(th(1))*sin(al(1)) a(1)*cos(th(1)); ...
      sin(th(1))   cos(th(1))*cos(al(1)) -cos(th(1))*sin(al(1)) a(1)*sin(th(1)); ...
               0              sin(al(1))             cos(al(1))            d(1); ...
               0                       0                      0              1];
A2 = [cos(th(2))  -sin(th(2))*cos(al(2))  sin(th(2))*sin(al(2)) a(2)*cos(th(2)); ...
      sin(th(2))   cos(th(2))*cos(al(2)) -cos(th(2))*sin(al(2)) a(2)*sin(th(2)); ...
               0              sin(al(2))             cos(al(2))            d(2); ...
               0                       0                      0              1];
A3 = [cos(th(3))  -sin(th(3))*cos(al(3))  sin(th(3))*sin(al(3)) a(3)*cos(th(3)); ...
      sin(th(3))   cos(th(3))*cos(al(3)) -cos(th(3))*sin(al(3)) a(3)*sin(th(3)); ...
               0              sin(al(3))             cos(al(3))            d(3); ...
               0                       0                      0              1];
A4 = [cos(th(4))  -sin(th(4))*cos(al(4))  sin(th(4))*sin(al(4)) a(4)*cos(th(4)); ...
      sin(th(4))   cos(th(4))*cos(al(4)) -cos(th(4))*sin(al(4)) a(4)*sin(th(4)); ...
               0              sin(al(4))             cos(al(4))            d(4); ...
               0                       0                      0              1];
A5 = [cos(th(5))  -sin(th(5))*cos(al(5))  sin(th(5))*sin(al(5)) a(5)*cos(th(5)); ...
      sin(th(5))   cos(th(5))*cos(al(5)) -cos(th(5))*sin(al(5)) a(5)*sin(th(5)); ...
               0              sin(al(5))             cos(al(5))            d(5); ...
               0                       0                      0              1];
A6 = [cos(th(6))  -sin(th(6))*cos(al(6))  sin(th(6))*sin(al(6)) a(6)*cos(th(6)); ...
      sin(th(6))   cos(th(6))*cos(al(6)) -cos(th(6))*sin(al(6)) a(6)*sin(th(6)); ...
               0              sin(al(6))             cos(al(6))            d(6); ...
               0                       0                      0              1];
           
A7 = [cos(th(7))  -sin(th(7))*cos(al(7))  sin(th(7))*sin(al(7)) a(7)*cos(th(7)); ...
      sin(th(7))   cos(th(7))*cos(al(7)) -cos(th(7))*sin(al(7)) a(7)*sin(th(7)); ...
               0              sin(al(7))             cos(al(7))            d(7); ...
               0                       0                      0              1];

           
g=[0 0 -9.80665 0];
Q= [0 -1 0 0
    1  0 0 0
    0  0 0 0
    0  0 0 0];
m1=5.70044;m2=3.22698;m3=4.31272;m4=2.07206;m5=2.24665;m6=1.60979;m7=0.54218;
x1=[-0.05117 0.07908 0.00086 1];%Column vector
x2=[0.00269 -0.00529 0.06845 1];
x3=[-0.07176 0.08149 0.00132 1];
x4=[0.00159 -0.01117 0.02618 1];
x5=[-0.01168 0.13111 0.00460 1];
x6=[0.00697  0.00600 0.06048 1];
x7=[0.005137 0.0009572 -0.06682 1];
%calculate gravity term
T10=A1;T20=A1*A2;T30=A1*A2*A3;T40=A1*A2*A3*A4;T50=A1*A2*A3*A4*A5;T60=A1*A2*A3*A4*A5*A6;T70=A1*A2*A3*A4*A5*A6*A7;
T21=A2;T31=A2*A3;T41=A2*A3*A4;T51=A2*A3*A4*A5;T61=A2*A3*A4*A5*A6;T71=A2*A3*A4*A5*A6*A7;
T32=A3;T42=A3*A4;T52=A3*A4*A5;T62=A3*A4*A5*A6;T72=A3*A4*A5*A6*A7;
T43=A4;T53=A4*A5;T63=A4*A5*A6;T73=A4*A5*A6*A7;
T54=A5;T64=A5*A6;T74=A5*A6*A7;
T65=A6;T75=A6*A7;
T76=A7;

u11=Q*T10;      u21=Q*T20;      u31=Q*T30;     u41=Q*T40;      u51=Q*T50;     u61=Q*T60;     u71=Q*T60;
u22=T10*Q*T21;  u32=T10*Q*T31;  u42=T10*Q*T41; u52=T10*Q*T51;  u62=T10*Q*T61; u72=T10*Q*T71;
u33=T20*Q*T32;  u43=T20*Q*T42;  u53=T20*Q*T52; u63=T20*Q*T62;  u73=T20*Q*T72;
u44=T30*Q*T43;  u54=T30*Q*T53;  u64=T30*Q*T63; u74=T30*Q*T73;
u55=T40*Q*T54;  u65=T40*Q*T64;  u75=T40*Q*T74;
u66=T50*Q*T65;  u76=T50*Q*T75;
u77=T60*Q*T76;

g1=-(m1*g*u11*x1'+m2*g*u21*x2'+m3*g*u31*x3'+m4*g*u41*x4'+m5*g*u51*x5'+m6*g*u61*x6'+m7*g*u71*x7');
g2=-(m2*g*u22*x2'+m3*g*u32*x3'+m4*g*u42*x4'+m5*g*u52*x5'+m6*g*u62*x6'+m7*g*u72*x7');
g3=-(m3*g*u33*x3'+m4*g*u43*x4'+m5*g*u53*x5'+m6*g*u63*x6'+m7*g*u73*x7');
g4=-(m4*g*u44*x4'+m5*g*u54*x5'+m6*g*u64*x6'+m7*g*u74*x7');
g5=-(m5*g*u55*x5'+m6*g*u65*x6'+m7*g*u75*x7');
g6=-(m6*g*u66*x6'+m7*g*u76*x7');
g7=-(m7*g*u77*x7');

%与加速度相关的惯量矩阵项
%calculate inertia term D(q),symmetric, bounded, positive definite
% D=[D11 D12 D13 D14 D15 D16 D17
%    D21 D22 D23 D24 D25 D26 D27
%    D31 D32 D33 D34 D35 D36 D37
%    D41 D42 D43 D44 D45 D46 D47 
%    D51 D52 D53 D54 D55 D56 D57
%    D61 D62 D63 D64 D65 D66 D67 
%    D71 D72 D73 D74 D75 D76 D77];


t1=g1;
t2=g2;
t3=g3;
t4=g4;
t5=g5;
t6=g6;
t7=g7;

G3=[t1;t2;t3;t4;t5;t6;t7];









