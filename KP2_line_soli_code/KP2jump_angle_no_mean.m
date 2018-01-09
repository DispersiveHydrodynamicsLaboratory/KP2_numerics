function [R1,R2,f,g] = KP2jump_angle_no_mean(R1m,R1p,R2m,R2p,m,choice)
% Calculates evolution of a jump in R1 and R2
% Where the jump occurs over a line with slope m
% Choice determines which solution branch is taken

%% Calculate slope of line perpendicular to m
mrec = -1/m;

%% Calculate characteristic speeds
W1 = @(R1) -m/(3*sqrt(1+m^2)) * (R1 + 2/m) * (R1 - 2*R2m);
W1m = W1(R1m);
W1p = W1(R1p);

W2 = @(R2) -m/(3*sqrt(1+m^2)) * (R2 - 2/m) * (R2 - 2*R1p);
W2m = W2(R2m);
W2p = W2(R2p);

%% Compute simple wave solutions
% if choice == 2
%     f = @(xi) R2m + mrec + sqrt( (R2m-mrec)^2 - 3*sqrt(1+m^2)/m * xi );
%     g = @(xi) R1p - mrec - sqrt( (R1p-mrec)^2 - 3*sqrt(1+m^2)/m * xi );
% elseif choice ==1
    f = @(xi) R2m + mrec - sqrt( (R2m-mrec)^2 - 3*sqrt(1+m^2)/m *xi ) +4.5;
    g = @(xi) R1p - mrec - sqrt( (R1p+mrec)^2 - 3*sqrt(1+m^2)/m *xi );
% end

%% Stitch together solutions for R1 and R2
if W1m>W1p
    warning('Reversed R1');
    c1m = W1p;
    c1p = W1m;
else
    c1m = W1m;
    c1p = W1p;
end
if W2m>W2p
    warning('Reversed R2');
    c2m = W2p;
    c2p = W2m;
else
    c2m = W2m;
    c2p = W2p;
end
alp = @(x,y) x + mrec*y;
R1 = @(x,y,t) R1m           .* ( alp(x,y) <= c1m*t ) + ...
              f(alp(x,y)/t) .* ( alp(x,y)  > c1m*t ) .* ( alp(x,y) <  c1p*t ) + ...
              R1p                                    .* ( alp(x,y) >= c1p*t );
R2 = @(x,y,t) R2m           .* ( alp(x,y) <= c2m*t ) + ...
              g(alp(x,y)/t) .* ( alp(x,y)  > c2m*t ) .* ( alp(x,y) <  c2p*t ) + ...
              R2p                                    .* ( alp(x,y) >= c2p*t );

%% Figure for debugging R1 & R2
figure(2); clf;
t = 100; x_tra = 1;
x = (t*(-x_tra+min([W1m W1p W2m W2p 0]))):0.001:(t*(x_tra+max([W1m W1p W2m W2p 0]))); 
sa = @(x,y,t) 1/2 * (R1(x,y,t) + R2(x,y,t));
a  = @(x,y,t) sa(x,y,t).^2;
q  = @(x,y,t) 1/2 * (R1(x,y,t) - R2(x,y,t));
    subplot(2,2,1:2);
        plot(alp(x,0)/t,R1(x,0,t),'b.')
        hold on;
        plot([W1m W1m],[R1m R1p],'b--',...
             [W1p W1p],[R1m R1p],'k--',...
             [W1m W1p],[R1m R1m],'b--',...
             [W1m W1p],[R1p R1p],'k--');
%          title('R_1')
%     subplot(2,2,2);
        plot(alp(x,0)/t,R2(x,0,t),'r.')
        hold on;
        plot([W2m W2m],[R2m R2p],'r-.',...
             [W2p W2p],[R2m R2p],'m-.',...
             [W2m W2p],[R2m R2m],'r-.',...
             [W2m W2p],[R2p R2p],'m-.');
        title('R_1 (blue) and R_2 (red)');
    subplot(2,2,3);
        plot(alp(x,0),a(x,0,t),'.');
        title('a');
    subplot(2,2,4);
        plot(alp(x,0),q(x,0,t),'.');
        title('q');
  disp('fiddle');        
          