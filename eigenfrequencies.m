% calculate the eigenfrequencies of the system based on the 2-node bar model
clear all;
clc;

syms omega;

L   = 12;        % inches, shaker-start-tail span
w   = 1;         % inch, width (all members)
h   = 1/8;       % inch, fuselage and wing thickness
h_E = 1/4;       % inch, elevator thickness
h_R = 0.040;     % inch, rudder thickness
E   = 10175000;  % lb/in^2 (70 GPA), modulus of elasticity
r   = 0.0002505; % lb * sec^2/in^4 (mass density)

A = w*h;
Izz = w*h^3/12;

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
          0       0     -642*L -73*L^2 -1458*L 172*L^2];
M2 = M2 + [0 0 0 0 0  0;
           0 0 0 0 0  0;
           0 0 0 0 0  0;
           0 0 0 0 0  0;
           0 0 0 0 MT ST;
           0 0 0 0 ST IT];


K2 = ck2*[24  6*L   -24  6*L    0   0;
          6*L 2*L^2 -6*L L^2    0   0;
         -24 -6*L    48  0     -24  6*L;
          6*L L^2    0   4*L^2 -6*L L^2;
          0   0     -24 -6*L    24 -6*L;
          0   0      6*L L^2   -6*L 2*L^2];

% first node is fixed, so it's removed
M2 = M2(3:end, 3:end);
K2 = K2(3:end, 3:end);

f = @(omega) det(K2 - (omega^2)*M2); % eigenfrequencies are where the determinant is zero

% scroll through all the frequencies and figure out where it's crosses zero
prev = 1;
count = 0;
freq = 0.01;
while count < 3
  this = f(freq);
  if prev * this < 1 % change sign
    fprintf('%f Hz\n', freq/(2*pi));
    count = count + 1;
  end
  prev = this;
  freq = freq + 0.005;
end
