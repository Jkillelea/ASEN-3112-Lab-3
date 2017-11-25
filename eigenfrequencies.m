% calculate the eigenfrequencies of the system based on the 2-node bar model
clear all;
clc;

L   = 12;        % inches, shaker-start-tail span
w   = 1;         % inch, width (all members)
h   = 1/8;       % inch, fuselage and wing thickness
h_E = 1/4;       % inch, elevator thickness
h_R = 0.040;     % inch, rudder thickness
E   = 10175000;  % lb/in^2 (70 GPA), modulus of elasticity
r   = 0.0002505; % lb * sec^2/in^4 (mass density)

A   = w*h;
Izz = w*h^3/12;

disp('============== 2 Element Model  ==============');
% idk what these represent
cm2 = r*A*L/100800;
ck2 = 4*E*Izz/L^3;

MT = 1.131*r;
ST = 0.5655*r;
IT = 23.124*r;

% mass and stiffness matricies
M2 = cm2*[19272  1458*L  5928  -642*L    0     0;
          1458*L 172*L^2 642*L -73*L^2   0     0;
          5928   642*L   38544  0       5928  -642*L;
          -642*L -73*L   0      344*L^2 642*L -73*L^3;
          0       0      5928   642*L   19272 -1458*L;
          0       0     -642*L -73*L^2 -1458*L 172*L^2] ...
...
       + [0       0      0      0       0      0;
          0       0      0      0       0      0;
          0       0      0      0       0      0;
          0       0      0      0       0      0;
          0       0      0      0       MT     ST;
          0       0      0      0       ST     IT];

K2 = ck2*[24  6*L   -24  6*L    0   0;
          6*L 2*L^2 -6*L L^2    0   0;
         -24 -6*L    48  0     -24  6*L;
          6*L L^2    0   4*L^2 -6*L L^2;
          0   0     -24 -6*L    24 -6*L;
          0   0      6*L L^2   -6*L 2*L^2];

% first node is fixed, so it's removed
M2 = M2(3:end, 3:end);
K2 = K2(3:end, 3:end);

[v, d] = eig(M2^-1 * K2);
% frequencies are the sqrt of the eigenvalues (sum because they're returned as diagonals of a matrix)
d = sum(sqrt(d));

for i = length(d):-1:(length(d) - 2) % lowest frequency is the last in the vector, so go in reverse
  fprintf('Frequency %d: %f Hz\n', (length(d) - i + 1), d(i)/(2*pi));
  fprintf('u: ');
  disp(v(:, i));
end


disp('============== 4 Element Model  ==============');
% 4 element model
cm4 = r*A*L/806400;
ck4 = 8*E*Izz/L^3;

M4 = cm4 * [77088   2916*L  23712   -1284*L  0        0        0         0       0       0;
            2916*L  172*L^2 1284*L  -73*L^2  0        0        0         0       0       0;
            23712   1284*L  154176  0        23712  -1284*L    0         0       0       0;
           -1284*L -73*L^2  0       344*L^2  1284*L  -73*L^2   0         0       0       0;
            0       0       23712   1284*L   154176   0        23712    -1284*L  0       0;
            0       0      -1284*L -73*L^2   0        344*L^2  1284*L   -73*L^2  0       0;
            0       0       0       0        23712    1284*L   154176    0       23712  -1284*L;
            0       0       0       0       -1284*L  -73*L^2   0         344*L^2 1284*L -73*L^2;
            0       0       0       0        0        0        23712     1284*L  77088  -2916*L;
            0       0       0       0        0        0       -1284*L   -73*L^2 -2916*L  172*L^2] ...
      + ...
     [0       0       0       0        0        0       0         0       0       0;
      0       0       0       0        0        0       0         0       0       0;
      0       0       0       0        0        0       0         0       0       0;
      0       0       0       0        0        0       0         0       0       0;
      0       0       0       0        0        0       0         0       0       0;
      0       0       0       0        0        0       0         0       0       0;
      0       0       0       0        0        0       0         0       0       0;
      0       0       0       0        0        0       0         0       0       0;
      0       0       0       0        0        0       0         0       MT     ST;
      0       0       0       0        0        0       0         0       ST     IT];


K4 = ck4*[96    12*L   -96   12*L   0    0      0     0      0      0;
          12*L  2*L^2  -12*L L^2    0    0      0     0      0      0;
         -96   -12*L    192   0    -96   12*L   0     0      0      0;
          12*L  L^2     0    4*L^2 -12*L L^2    0     0      0      0;
          0     0      -96  -12*L   192  0     -96    12*L   0      0;
          0     0       12*L L^2    0    4*L^2 -12*L  L^2    0      0;
          0     0       0    0     -96  -12*L   192   0     -96     12*L;
          0     0       0    0      12*L L^2    0     4*L^2 -12*L   L^2;
          0     0       0    0      0    0     -96   -12*L   96    -12*L;
          0     0       0    0      0    0      12*L  L^2   -12*L   2*L^2];

M4 = M4(3:end, 3:end);
K4 = K4(3:end, 3:end);

[v, d] = eig(M4^-1 * K4);
% frequencies are the sqrt of the eigenvalues (sum because they're returned as diagonals of a matrix)
d = sqrt(sum(d))';

for i = length(d):-1:(length(d) - 2) % lowest frequency is the last in the vector, so go in reverse
  fprintf('Frequency %d: %f Hz\n', (length(d) - i + 1), d(i)/(2*pi));
  fprintf('u: ');
  disp(v(:, i));
end
