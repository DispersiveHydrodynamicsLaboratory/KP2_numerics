disp('Running hello...');
q = strsplit(pwd,filesep);
gpath = [strjoin(q(1:end-1),filesep),filesep,'solver_GPU'];
addpath(gpath);
matlab_check;
save_dir = '/rc_scratch/mima6446/';
A = magic(3);
save([save_dir,'test.mat'],'A');
