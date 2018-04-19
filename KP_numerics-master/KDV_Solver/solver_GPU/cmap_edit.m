function [ cmap2 ] = cmap_edit( cred, cmin, cmax, cmap )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
ctop = cmap(end,:);
cmid = cmap(256/2+1,:);
cbot = cmap(1,:);
cper = round((cred-cmin)/(cmax-cmin)*100)/100; % Cred as a fraction of clims
csplit = round(300*cper);
cmap2 = zeros(300,3);

    cmap2(1:csplit,1) = interp1(linspace(1,csplit,129),cmap(1:129,1),1:csplit);
    cmap2(1:csplit,2) = interp1(linspace(1,csplit,129),cmap(1:129,2),1:csplit);
    cmap2(1:csplit,3) = interp1(linspace(1,csplit,129),cmap(1:129,3),1:csplit);

    cmap2(csplit+1:end,1) = interp1(linspace(1,300-csplit,129),cmap(129:end,1),1:(300-csplit));
    cmap2(csplit+1:end,2) = interp1(linspace(1,300-csplit,129),cmap(129:end,2),1:(300-csplit));
    cmap2(csplit+1:end,3) = interp1(linspace(1,300-csplit,129),cmap(129:end,3),1:(300-csplit));


end