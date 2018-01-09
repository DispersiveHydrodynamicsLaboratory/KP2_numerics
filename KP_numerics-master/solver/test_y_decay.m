%% Script for playing with different y-reducers

%% Want a function that:
%   is 1 for [0, +Ly-b]
%   goes from 1 to 0 from [Ly-b, Ly-c]
%   is 0 for [Ly-c, Ly]
Ly = 20;
b = 10;
c = 5;

% Tukey window (thanks Wikipedia!)
w = @(Y) repmat(tukeywin(size(Y,1),(b-c)/Ly),[1,size(Y,2)]);

% %% Even reflection of function over y-axis
% w = @(Y) (w0(Y) + w0(-Y))/2;

% Test 
y = -Ly:0.1:Ly;
Y = repmat(y',[1 5]);
figure(3); clf;
    plot(y,w(Y));