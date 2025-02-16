% dzhanibekov effect simulation
clear; clc; close all;

% moments of inertia (i1 < i2 < i3; i2 is unstable)
I1 = 1; I2 = 2; I3 = 3;

% initial conditions: rotate about the unstable (i2) axis
omega0 = [0; 1; 0];        % angular velocity in body frame
q0     = [1; 0; 0; 0];     % unit quaternion (no rotation)
x0     = [omega0; q0];

tspan = [0 20];
opts  = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t,x] = ode45(@(t,x) rigidODE(t,x,I1,I2,I3), tspan, x0, opts);

% setup figure and cube vertices
figure; axis equal; grid on; view(3);
xlabel('x'); ylabel('y'); zlabel('z');
xlim([-2 2]); ylim([-2 2]); zlim([-2 2]);
cube = [ -0.5 -0.5 -0.5;
          0.5 -0.5 -0.5;
          0.5  0.5 -0.5;
         -0.5  0.5 -0.5;
         -0.5 -0.5  0.5;
          0.5 -0.5  0.5;
          0.5  0.5  0.5;
         -0.5  0.5  0.5];
faces = [1 2 3 4;
         5 6 7 8;
         1 2 6 5;
         2 3 7 6;
         3 4 8 7;
         4 1 5 8];

% animate rotation using computed quaternions
for k = 1:length(t)
    q = x(k,4:7)';                   % current orientation
    R = quat2rotm(q');                % get rotation matrix
    rotatedCube = (R * cube')';       % rotate vertices
    cla; patch('Vertices', rotatedCube, 'Faces', faces, ...
              'FaceColor','cyan','FaceAlpha',0.8);
    title(sprintf('t = %.2f sec', t(k)));
    drawnow; pause(0.02);
end

% --- helper functions ---

function dx = rigidODE(~, x, I1, I2, I3)
    omega = x(1:3);
    q     = x(4:7);
    % torque-free euler equations
    domega = [((I2-I3)/I1)*omega(2)*omega(3);
              ((I3-I1)/I2)*omega(3)*omega(1);
              ((I1-I2)/I3)*omega(1)*omega(2)];
    % quaternion derivative (q = [q0,q1,q2,q3], q0 is scalar)
    qdot = 0.5 * [ -q(2)*omega(1) - q(3)*omega(2) - q(4)*omega(3);
                    q(1)*omega(1) - q(4)*omega(2) + q(3)*omega(3);
                    q(4)*omega(1) + q(1)*omega(2) - q(2)*omega(3);
                   -q(3)*omega(1) + q(2)*omega(2) + q(1)*omega(3)];
    dx = [domega; qdot];
end

function R = quat2rotm(q)
    % converts quaternion [q0,q1,q2,q3] to rotation matrix
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    R = [1-2*(q2^2+q3^2), 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2);
         2*(q1*q2+q0*q3), 1-2*(q1^2+q3^2), 2*(q2*q3-q0*q1);
         2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 1-2*(q1^2+q2^2)];
end
