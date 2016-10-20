close all
clear -regexp
clc

tic           % start timer

dt = 1 / 100; % sample time
r = 4;        % derivative to minimize in cost function (snap)
n = 2*r -1;   % order of desired trajectory
m = 3;        % number of pieces in trajectory (m+1 waypoints) 
d = 3;        % dimensions (d=3 means 3D)

save_flag = false;       % true to save output for plots
nc = 20;                 % number constraints between two keyframe
corridor_radious = 0.03; % radius of the constrain tube
corridor_start = 2;      % specify on wich waypoint you want to start the 
                         % corridor constraints
corridor_end = 4;        % specify on wich waypoint you want to end the                  
                         % corridor constraints (DOESN'T WORK FOR NOW)

% specify desired arrival times at keyframes (in seconds)
tDes = [0;    % keep the first equal to 0
        10; 
        20; 
        30]; 

% For each waypoints you need to create an object of m+1 columns and r+1 
% rows. In the first row write the desire position at each waypoints, while 
% in the next rows write the desire velocity, acceleration, jerk and snap.
% If 255 means that that condition is free (no constrained). Repeat for
% each dimension.
    
% X
posDes(:, :, 1) = [0 1.5 2.1 4.2; ...
                   0 255 255 0;   ... % initial and final velocity is 0
                   0 255 255 0;   ... % initial and final acceleration is 0
                   0 255 255 0;   ... % initial and final jerk is 0
                   0 255 255 0];      % initial and final snap is 0
% Y               
posDes(:, :, 2) = [0 0   0.9 1.2; ...
                   0 255 255 0;   ...
                   0 255 255 0;   ...
                   0 255 255 0;   ...
                   0 255 255 0];   

% Z
posDes(:, :, 3) = [0 0.3 0.4 1; ...
                   0 255 255 0; ...
                   0 255 255 0; ...
                   0 255 255 0; ...
                   0 255 255 0];
 
 
%% find trajectories for each dimension, nondimensionalized in time

xT = zeros(n+1, m, d); 
posDes_opt = zeros(r, m+1, d); 
 
for i=1 : d
   xT(:, :, i) = findTraj(r, n, m, i, tDes, posDes);
end


%% BETA corridor constraints 

ineqConst.nc = nc;
ineqConst.start = corridor_start * ones(ineqConst.nc, 1); 
ineqConst.delta = corridor_radious; 
ineqConst.dim = [1*ones(nc, 1), 2*ones(nc, 1), 3*ones(nc, 1)];
xTcorr = findTrajCorr(r, n, m, d, tDes, posDes, ineqConst);

%% evaluate QP traj

t = 0 : dt : tDes(m+1); % construct t vector 
t = t';
numDer = r;
pos = zeros(numDer+1, d, length(t));
posCorr = zeros(numDer+1, d, length(t));

% evaluate the piece-wise polynominal at each point in time
for i=1 : length(t)  
    [~, derivativesX] = evaluateTraj(t(i, 1), n, m, d, xT, tDes, ...
        numDer, []);
    [~, derivativesXCorr] = evaluateTraj(t(i, 1), n, m, d, xTcorr, tDes, ...
        numDer, []);
    
    pos(:, :, i) = evaluateTraj(t(i, 1), n, m, d, xT, tDes, numDer, ...
        derivativesX);
    posCorr(:, :, i) = evaluateTraj(t(i, 1), n, m, d, xTcorr, tDes, ...
        numDer, derivativesXCorr);
end

time = toc    % end timer


%% plot QP traj

% save data for plot
posx(:, 1) = pos(1, 1, :);
posy(:, 1) = pos(1, 2, :);
posz(:, 1) = pos(1, 3, :);
velx(:, 1) = pos(2, 1, :);
vely(:, 1) = pos(2, 2, :);
velz(:, 1) = pos(2, 3, :);
accx(:, 1) = pos(3, 1, :);
accy(:, 1) = pos(3, 2, :);
accz(:, 1) = pos(3, 3, :);
jerx(:, 1) = pos(4, 1, :);
jery(:, 1) = pos(4, 2, :);
jerz(:, 1) = pos(4, 3, :);
snax(:, 1) = pos(5, 1, :);
snay(:, 1) = pos(5, 2, :);
snaz(:, 1) = pos(5, 3, :);

posCorrx(:, 1) = posCorr(1, 1, :);
posCorry(:, 1) = posCorr(1, 2, :);
posCorrz(:, 1) = posCorr(1, 3, :);

figure
hold on
plot3(posx, posy, posz)
plot3(posCorrx, posCorry, posCorrz)
plot3(posDes(1, :, 1), posDes(1, :, 2), posDes(1, :, 3))
plot3(posDes(1, :, 1), posDes(1, :, 2), posDes(1, :, 3), 'k^')
hold off
legend('polynomial trajectory', 'corridor', 'trivial trajectory', 'setpoints')
title('position')

% uncomment to print velocity and acceleration
% figure
% hold on
% plot3(velx, vely, velz)
% hold off
% title('velocity')
% 
% figure
% hold on
% plot3(accx, accy, accz)
% hold off
% title('acceleration')


%% save DAT file for LaTeX
if(~save_flag)
    return;
end

% save .mat files
save('posx.mat', 'posx')
save('posy.mat', 'posy')
save('posz.mat', 'posz')
save('velx.mat', 'velx')
save('vely.mat', 'vely')
save('velz.mat', 'velz')
save('accx.mat', 'accx')
save('accy.mat', 'accy')
save('accz.mat', 'accz')

traj = zeros(length(posx), 2);
traj_x = zeros(length(posx), 2);
traj_z = zeros(length(posx), 2);
traj_xd = zeros(length(velx), 2);
traj_xdd = zeros(length(accx), 2);
traj_xddd = zeros(length(jerx), 2);
traj_xdddd = zeros(length(snax), 2);
for i=1 : length(traj)
    traj(i, 1) = posx(i);
    traj(i, 2) = posy(i);
    
    traj_x(i, 1) = t(i);    
    traj_x(i, 2) = posx(i);
    
    traj_z(i, 1) = t(i);    
    traj_z(i, 2) = posz(i);
    
    traj_xd(i, 1) = t(i);
    traj_xd(i, 2) = velx(i);
    
    traj_xdd(i, 1) = t(i);
    traj_xdd(i, 2) = accx(i); 
    
    traj_xddd(i, 1) = t(i);
    traj_xddd(i, 2) = jerx(i); 
    
    traj_xdddd(i, 1) = t(i);
    traj_xdddd(i, 2) = snax(i);     
end

traj_corr = zeros(length(posCorrx), 2);
for i=1 : length(traj_corr)
    traj_corr(i, 1) = posCorrx(i);
    traj_corr(i, 2) = posCorry(i);    
end

save('traj.dat', 'traj', '-ascii');
save('traj_x.dat', 'traj_x', '-ascii');
save('traj_z.dat', 'traj_z', '-ascii');
save('vel_x.dat', 'traj_xd', '-ascii');
save('acc_x.dat', 'traj_xdd', '-ascii');
save('jerk_x.dat', 'traj_xddd', '-ascii');
save('snap_x.dat', 'traj_xdddd', '-ascii');
save('traj_corr.dat', 'traj_corr', '-ascii');


desx = posDes(1, :, 1);
desy = posDes(1, :, 2);
desz = posDes(1, :, 3);
des_pos = [desx' desy'];
save('des_pos.dat', 'des_pos', '-ascii');
des_pos_z = [tDes desz'];
save('des_pos_z.dat', 'des_pos_z', '-ascii');
save('des_pos.dat', 'des_pos', '-ascii');
