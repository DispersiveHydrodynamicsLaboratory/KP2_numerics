function[savefile] = reformat_numerics_files(data_dir);
    loadfile = sprintf('%sparameters.mat',data_dir);
    load(loadfile);
    %% Calculating maximum time increment in saved data files...
    for ii=1:length(t)+1
        [fid,foo] = fopen(strcat(data_dir,num2str(ii,'%05d.mat')),'r');
        if fid == -1 % File does not exist
            tm = ii-1;
            disp(['Maximum time = ',num2str(t(tm))]);
            break;
        end
        fclose(fid);
    end
    if ii == length(t)
        disp(['Maximum time = ',num2str(t(tm))]);
    end
    % Get rid of larger t values
    t = t(1:tm);
    %% Format data structure for area matrix
    area_mat = zeros(length(t),Nz);
    
    %% Generate area matrix
	load(strcat(data_dir,num2str(0,'%05d')),'A_init');
    area_mat(1,:) = A_init;
    disp('Generating matrix of data...');
    for ti = 2:tm
        load(strcat(data_dir,num2str(ti,'%05d')),'A');
        area_mat(ti,:) = A;
    end

    %% Generate diameter, left_edges and right_edges matrices
    diam        = 2*sqrt(area_mat/pi);
    left_edges  = -diam/2;
    right_edges =  diam/2;
    disp('Saving...');
    savefile = [data_dir, 'conduit_edges.mat'];
    save(savefile,'area_mat','diam','left_edges','right_edges');
    disp(['Process has completed. File is saved to ', savefile]);