%{
Brian Knisely
ME523, HW1, P3
January 21, 2018

The purpose of this code is to compare five differencing schemes for
taking the first derivative with respect to x of cos(x) at x = 0.3
and calculating the corresponding errors. The schemes to be compared are
BDS1, CDS2, FDS3, BDS3, and CDS4.
%}

clear; close all; format compact; home;

n = 41;  % Number of iterations to compute
x = 0.3;  % Location to evaluate derivative
h = zeros(n, 1);  % Initialize array of spacings
h(1) = 2;  % Initialize first value of spacing to be 2 
for i = 1:n-1
    h(i+1) = h(i)/2;  % Reduce size of spacing by 50% each iteration
end


% Initialize arrays for each scheme
% The first column in each array holds its value, and the second column
% stores its error
bds1 = zeros(n, 2);  % 1st order backward difference
cds2 = zeros(n, 2);  % 2nd order central difference
fds3 = zeros(n, 2);  % 3rd order forward difference
bds3 = zeros(n, 2);  % 3rd order backward difference
cds4 = zeros(n, 2);  % 4th order central difference


for i = 1:n
    exact = -sin(0.3);  % Compute analytical derivative
    % Approximate derivative with each scheme
    fds1(i, 1) = (cos(x+h(i))-cos(x)) / h(i);
    bds1(i, 1) = (cos(x)-cos(x-h(i))) / h(i);
    cds2(i, 1) = (cos(x+h(i))-cos(x-h(i))) / (2*h(i));
    fds3(i, 1) = (-cos(x+2*h(i))+6*cos(x+h(i))-3*cos(x)-2*cos(x-h(i))) /...
        (6*h(i));
    bds3(i, 1) = (2*cos(x+h(i))+3*cos(x)-6*cos(x-h(i))+cos(x-2*h(i))) / ...
        (6*h(i));
    cds4(i, 1) = (-cos(x+2*h(i))+8*cos(x+h(i))-8*cos(x-h(i))+...
        cos(x-2*h(i))) / (12*h(i));
    % Compute errors for each scheme
    fds1(i, 2) = fds1(i, 1) - exact;
    bds1(i, 2) = bds1(i, 1) - exact;
    cds2(i, 2) = cds2(i, 1) - exact;
    fds3(i, 2) = fds3(i, 1) - exact;
    bds3(i, 2) = bds3(i, 1) - exact;
    cds4(i, 2) = cds4(i, 1) - exact;
end

figure(1);  % Create figure to show error
loglog(h, abs(fds1(:, 2)), 'r+-', h, abs(bds1(:, 2)), 'ko--', ...
    h, abs(cds2(:, 2)), 'rx:', h, abs(fds3(:, 2)), 'b*-', ...
    h, abs(bds3(:, 2)), 'ms--', h, abs(cds4(:, 2)), 'kp:'); 
xlabel('Spacing \it{h}'); ylabel('|Error|'); 
legend('FDS1', 'BDS1', 'CDS2', 'FDS3', 'BDS3', 'CDS4', ...
    'location', 'southwest'); 
set(gca, 'fontsize', 16);
xlim([1e-13, 10]);
set(gcf, 'outerposition', [50 50 850 650])

fid = fopen('results.txt', 'w');  % Open a file to store results

% Print error results
fprintf(fid,'%-2s|%-12s|%-12s|%-12s|%-12s|%-12s|%-12s|%-12s\n', ...
    'i', '      h', ' Error:FDS1', ' Error:BDS1', ' Error:CDS2', ...
    ' Error:FDS3', ' Error:BDS3', ' Error:CDS4');
for i = 1:n
    fprintf(fid,'%2.0f,%12.5e,%12.5e,%12.5e,%12.5e,%12.5e,%12.5e,%12.5e', ...
    i, h(i), fds1(i, 2), bds1(i, 2), cds2(i, 2), fds3(i, 2), bds3(i, 2), cds4(i, 2));
    fprintf(fid,'\n');
end

fclose(fid);  % Close the file