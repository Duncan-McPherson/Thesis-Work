function [phi, y, t, yd] = loadMunich(i)
    %% Load in Data from the Munich trials depending on Selection
    %Load Data From Files
    T_0 = 'Munich Data\Feedstep150_N0_H55.sco';
    T_con = 'Munich Data\Feedstep150_N300_H55.sco';
    T_sweep = 'Munich Data\Feedstep150_Nsweep300_H55.sco';
    
    tmp1 = load([T_con(1:end-3) 'mat']);
    t_con = tmp1.t_in;
    tmp1 = tmp1.AI;
    con_shaker = tmp1(:,1)*1000;
    [Fq_con_shak, Gy_con_shak] = fftspec(con_shaker', 500, 4e-4);
    
    tmp1 = load([T_sweep(1:end-3) 'mat']);
    t_sweep = tmp1.t_in;
    tmp1 = tmp1.AI;
    sweep_shaker = tmp1(:,1)*1000;
    [Fq_sweep_shak, Gy_sweep_shak] = fftspec(sweep_shaker', 500, 4e-4);
    
    [T_0_data] = read_tncscope_file(T_0);
    [T_con_data] = read_tncscope_file(T_con);
    [T_sweep_data] = read_tncscope_file(T_sweep);
    
    %Set Inital Range Around The Trajectory
    Ns = 1880;
    Ne = 4600;
    
    %Collect And Store Data From Inital Trajectory
    v_r = T_0_data.channel(2).values;
    v_m = T_0_data.channel(1).values;
    v_r_0 = v_r(Ns:Ne)*0.001/60; %% mm/min --> m/sec
    v_m_0 = v_m(Ns:Ne)*0.001/60; %% mm/min --> m/sec
    
    %Set Time Across The Range
    dt = 0.003;
    t = ([0:(Ne-Ns)]*dt)';
    
    %Select The Same Range In The Other Datasets
    v_r1 = T_con_data.channel(2).values;
    [~,offset_con] = max(xcorr(v_r,v_r1));
    offset_con = length(v_r1) - offset_con;
    v_m1 = T_con_data.channel(1).values;
    v_r_con = v_r1(Ns+offset_con:Ne+offset_con)*0.001/60; %% mm/min --> m/sec
    v_m_con = v_m1(Ns+offset_con:Ne+offset_con)*0.001/60; %% mm/min --> m/sec
    
    v_r2 = T_sweep_data.channel(2).values;
    [~,offset_sweep] = max(xcorr(v_r,v_r2));
    offset_sweep = length(v_r2) - offset_sweep;
    v_m2 = T_sweep_data.channel(1).values;
    v_r_sweep = v_r2(Ns+offset_sweep:Ne+offset_sweep)*0.001/60; %% mm/min --> m/sec
    v_m_sweep = v_m2(Ns+offset_sweep:Ne+offset_sweep)*0.001/60; %% mm/min --> m/sec
    
    %Create data vecotors as per the method setup
    deadBand = 0.0001; %Friction deadband

    switch i
        case 1
            y = v_m_0(2:end);
            phi = [v_m_0(1:end-1), v_r_0(1:end-1), -1*pv(v_m_0(1:end-1), deadBand), -1*nv(v_m_0(1:end-1), deadBand)];
            t = t(1:end-1);
            yd = zeros(1,length(y));
        case 4
            y = v_m_con(2:end);
            phi = [v_m_con(1:end-1), v_r_con(1:end-1), -1*pv(v_m_con(1:end-1), deadBand), -1*nv(v_m_con(1:end-1), deadBand)];
            t = t(1:end-1);
            yd = con_shaker;
        case 13
            y = v_m_sweep(2:end);
            phi = [v_m_sweep(1:end-1), v_r_sweep(1:end-1), -1*pv(v_m_sweep(1:end-1), deadBand), -1*nv(v_m_sweep(1:end-1), deadBand)];
            t = t(1:end-1);
            yd = sweep_shaker;
        otherwise
            warning("No data selected")
            y = 0; phi = 0; t = 0;
    end
end