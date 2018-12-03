%%%%%%%%%%%%%%%%% ANALYSING EIGENVALUES OF H %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

% last H was saved to file in 'basic_monodomain_1D', need to load it 1st
load('H.mat');  % also loads variables: Nx, Nt, N_tot

% assume we have 2 separate lists of space and time indices
% make other two of the same length ... 

% assume we are not choosing from the 1st & last time step -> all values
% are double, want to incorporate them both to get connection in time

zero_space_ind_1 = 16:20;
zero_time_ind_1 = 80:83;

zero_space_ind_2 = 13:17;
zero_time_ind_2 = 20:23;

one_space_ind = 4:8;
one_time_ind = 105:108;

% eigenvalues will come in sets, # of points in time: # points in set
% # of points in space: # of sets
no_of_sets = length(one_space_ind);
set_size = length(one_time_ind);

l_space = length(one_space_ind);
l_time = length(one_time_ind);

wave_space_ind_1 = 10:14;
wave_time_ind_1 = 77:80;

wave_space_ind_2 = 8:12;
wave_time_ind_2 = 30:33;

% translate time indices according to DG
zero_time_ind_1_DG = zeros(1,2*l_time);
zero_time_ind_2_DG = zeros(1,2*l_time);


one_time_ind_DG = zeros(1,2*l_time);
wave_1_time_ind_DG = zeros(1,2*l_time);
wave_2_time_ind_DG = zeros(1,2*l_time);



count = 1;
for k=1:l_time
    zero_time_ind_1_DG(count:count+1) = 2*zero_time_ind_1(k)-3:2*zero_time_ind_1(k)-2;
    zero_time_ind_2_DG(count:count+1) = 2*zero_time_ind_2(k)-3:2*zero_time_ind_2(k)-2;

    one_time_ind_DG(count:count+1) = 2*one_time_ind(k)-3:2*one_time_ind(k)-2;    
    wave_time_ind_1_DG(count:count+1) = 2*wave_time_ind_1(k)-3:2*wave_time_ind_1(k)-2;    
    wave_time_ind_2_DG(count:count+1) = 2*wave_time_ind_2(k)-3:2*wave_time_ind_2(k)-2; 
  count = count+2;
end

% draw points in graph of u ... so that we can see
figure; hold on;
mesh(x, t, u_mat'); hold on;
title('points chosen for submatrices');
xlabel('space');
ylabel('time');
for i=1:l_space
    for j=1:l_time
        plot3(x(zero_space_ind_1(i)),t(zero_time_ind_1(j)),u_mat(zero_space_ind_1(i),zero_time_ind_1(j)),'.r','markersize',10);
        plot3(x(zero_space_ind_2(i)),t(zero_time_ind_2(j)),u_mat(zero_space_ind_2(i),zero_time_ind_2(j)),'.r','markersize',10);
        
        plot3(x(one_space_ind(i)),t(one_time_ind(j)),u_mat(one_space_ind(i),one_time_ind(j)),'.r','markersize',10);
        plot3(x(wave_space_ind_1(i)),t(wave_time_ind_1(j)),u_mat(wave_space_ind_1(i),wave_time_ind_1(j)),'.r','markersize',10);
        plot3(x(wave_space_ind_2(i)),t(wave_time_ind_2(j)),u_mat(wave_space_ind_2(i),wave_time_ind_2(j)),'.r','markersize',10);

        %mesh(x(i), t(j), u_mat(i,j));
        hold on; 
    end
end

zero_sp_t_ind_1 = zeros(1, 2*l_space*l_time);
zero_sp_t_ind_2 = zeros(1, 2*l_space*l_time);

one_sp_t_ind = zeros(1, 2*l_space*l_time);
wave_sp_t_ind_1 = zeros(1, 2*l_space*l_time);
wave_sp_t_ind_2 = zeros(1, 2*l_space*l_time);

count = 1;

for j=1:2*l_time
    for i=1:l_space
        zero_sp_t_ind_1(count) = zero_time_ind_1_DG(j)*Nx + zero_space_ind_1(i);
        zero_sp_t_ind_2(count) = zero_time_ind_2_DG(j)*Nx + zero_space_ind_2(i);

        one_sp_t_ind(count) = one_time_ind_DG(j)*Nx + one_space_ind(i);
        wave_sp_t_ind_1(count) = wave_time_ind_1_DG(j)*Nx + wave_space_ind_1(i);
        wave_sp_t_ind_2(count) = wave_time_ind_2_DG(j)*Nx + wave_space_ind_2(i);

        count = count+1;
    end
end


H_zero_1 = H(zero_sp_t_ind_1, zero_sp_t_ind_1);
H_zero_2 = H(zero_sp_t_ind_2, zero_sp_t_ind_2);


H_one = H(one_sp_t_ind,one_sp_t_ind);
H_wave_1 = H(wave_sp_t_ind_1,wave_sp_t_ind_1);
H_wave_2 = H(wave_sp_t_ind_2,wave_sp_t_ind_2);

% u_sub_zero = u_mat(zero_space_ind, zero_time_ind);
% u_sub_one = u_mat(one_space_ind, one_time_ind);
% u_sub_wave_1 = u_mat(wave_space_ind_1, wave_time_ind_1);
% u_sub_wave_2 = u_mat(wave_space_ind_2, wave_time_ind_2);

[V_zero_1, D_zero_1] = eig(full(H_zero_1));
[V_zero_2, D_zero_2] = eig(full(H_zero_2));


[V_wave_1, D_wave_1] = eig(full(H_wave_1));
[V_wave_2, D_wave_2] = eig(full(H_wave_2));
[V_one, D_one] = eig(full(H_one));

%%%%%%%%%%%%%%%%%%%%% PLOT EIGENVECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_sm = linspace(0,l_space-1, l_space);
y_sm = linspace(0, l_time-1, 2*l_time);

%%%%%%%%%%%%%%%%% PLOT EIGENVALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pts = 1:l_space*2*l_time;
% 
figure;
plot(diag(D_zero_1), 'o'); hold on;
plot(diag(D_zero_2), 'o'); hold on;
plot(diag(D_one), 'o'); hold on;
plot(diag(D_wave_1), '*'); hold on;
plot(diag(D_wave_2), '*');
legend('ev zero zone 1', 'ev zero zone 2', 'ev one zone', 'ev wave front 1', 'ev wave front 2');
title('eigenvalues');

norm_diff = norm(D_zero_1 - D_zero_2);
%%%%% Grouping Eigenvectors 

tol = 10^(-17);
diff_norm_V_zero = zeros(size(V_zero_1));
diff_norm_ind_V_zero = zeros(size(V_zero_1));

pair_count = 0;

for it=1:size(V_zero_1,1)
    for it_2=it+1:size(V_zero_1,1)
        %temp_1 = norm(V_wave_2(:,it) - V_zero_1(:,it_2));
        %temp_2 = norm(V_wave_2(:,it) + V_zero_1(:,it_2));
        
        temp_com_conj = norm(real(V_wave_2(:,it)) - real(V_wave_2(:,it_2)));
        diff_norm_V_zero(it, it_2) = min(temp_1, temp_2);
        % if(diff_norm_V_zero(it, it_2) < tol)
        if(temp_com_conj < tol)
            pair_count = pair_count + 1;
            % diff_norm_ind_V_zero(it, it_2) = 1;
            break;
        end
    end
    
end

% pairs of eigenvectors that seem to be very similar
sparse(diff_norm_ind_V_zero);
nnz(diff_norm_ind_V_zero)

%%%%%%%%%%%% SORTED PLOTS BY EIGENVALUES %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% eigenfunctions are complex, for now plot them on top of each other
% where fact that one is imaginary part is indicated by colour + linstyle
% for k=1:2   %l_space
%     it=1;
%     
%     fig = figure;
%     suptitle(['one zone ', num2str(k), 'th set, REAL PART']);
%     for l=(k-1)*l_space*l_time+1:k*l_time*l_space
%         Evec = -1*ones(l_space, 2*l_time);
%         for ts=1:2*l_time
%             Evec(:, ts) = V_one(l_space*(ts-1)+1:l_space*ts,l);
%         end
%         subplot(4, ceil(l_space*l_time/4), it); 
%         mesh(x_sm, y_sm, real(Evec)', 'EdgeColor','b'); hold on;
%         %mesh(x_sm, y_sm, imag(Evec)', 'EdgeColor','g');
%         %legend('real part', 'imaginary part');
%         title(['EigVal ' num2str(D_one(l,l))]);
%         xlabel('space');
%         ylabel('time');
%         it = it+1;
%     end
%     %filename = ['eigvec_one_zone_set_', num2str(k), '_real_part'];
%     %saveas(fig, filename);
% end

% for k=1:2   %l_space
%     it=1;
%     
%     fig = figure;
%     suptitle(['wave 1 area ', num2str(k), 'th set, REAL PART']);
%     for l=(k-1)*l_space*l_time+1:k*l_time*l_space
%         Evec = -1*ones(l_space, 2*l_time);
%         for ts=1:2*l_time
%             Evec(:, ts) = V_wave_1(l_space*(ts-1)+1:l_space*ts,l);
%         end
%         subplot(4, ceil(l_space*l_time/4), it); 
%         mesh(x_sm, y_sm, real(Evec)', 'EdgeColor','b'); hold on;
%         %mesh(x_sm, y_sm, imag(Evec)', 'EdgeColor','g');
%         %legend('real part', 'imaginary part');
%         title(['EigVal ' num2str(D_wave_1(l,l))]);
%         xlabel('space');
%         ylabel('time');
%         it = it+1;
%     end
%     %filename = ['eigvec_wave_1_zone_set_', num2str(k), '_real_part'];
%     %saveas(fig, filename);
% end
% 
% for k=1:2   %l_space
%     it=1;
%     
%     fig = figure;
%     suptitle(['wave 2 area ', num2str(k), 'th set, REAL PART']);
%     for l=(k-1)*l_space*l_time+1:k*l_time*l_space
%         Evec = -1*ones(l_space, 2*l_time);
%         for ts=1:2*l_time
%             Evec(:, ts) = V_wave_2(l_space*(ts-1)+1:l_space*ts,l);
%         end
%         subplot(4, ceil(l_space*l_time/4), it); 
%         mesh(x_sm, y_sm, real(Evec)', 'EdgeColor','b'); hold on;
%         %mesh(x_sm, y_sm, imag(Evec)', 'EdgeColor','g');
%         %legend('real part', 'imaginary part');
%         title(['EigVal ' num2str(D_wave_2(l,l))]);
%         xlabel('space');
%         ylabel('time');
%         it = it+1;
%     end
%     %filename = ['eigvec_wave_2_zone_set_', num2str(k), '_real_part'];
%     %saveas(fig, filename);
% end
% 

for k=1:2
    it=1;
    
    fig = figure;
    suptitle(['zero zone 1 ', num2str(k), 'th set, REAL PART']);
    for l=(k-1)*l_space*l_time+1:k*l_time*l_space
        Evec = -1*ones(l_space, 2*l_time);
        for ts=1:2*l_time
            Evec(:, ts) = V_zero_1(l_space*(ts-1)+1:l_space*ts,l);
        end
        subplot(4, ceil(l_space*l_time/4), it); 
        mesh(x_sm, y_sm, real(Evec)', 'EdgeColor','b'); hold on;
        %mesh(x_sm, y_sm, imag(Evec)', 'EdgeColor','g');
        %legend('real part', 'imaginary part');
        title(['EigVal ' num2str(D_zero_1(l,l))]);
        xlabel('space');
        ylabel('time');
        it = it+1;
    end
    %filename = ['eigvec_zero_zone_set_', num2str(k), '_real_part'];
    %saveas(fig, filename);
end

for k=1:2
    it=1;
    
    fig = figure;
    suptitle(['zero zone 2 ', num2str(k), 'th set, REAL PART']);
    for l=(k-1)*l_space*l_time+1:k*l_time*l_space
        Evec = -1*ones(l_space, 2*l_time);
        for ts=1:2*l_time
            Evec(:, ts) = V_zero_2(l_space*(ts-1)+1:l_space*ts,l);
        end
        subplot(4, ceil(l_space*l_time/4), it); 
        mesh(x_sm, y_sm, real(Evec)', 'EdgeColor','b'); hold on;
        %mesh(x_sm, y_sm, imag(Evec)', 'EdgeColor','g');
        %legend('real part', 'imaginary part');
        title(['EigVal ' num2str(D_zero_2(l,l))]);
        xlabel('space');
        ylabel('time');
        it = it+1;
    end
    %filename = ['eigvec_zero_zone_set_', num2str(k), '_real_part'];
    %saveas(fig, filename);
end
% 
% for k=1:l_space
%     it=1;
%     
%     fig = figure;
%     suptitle(['Wave Front 1 ', num2str(k), 'th set']);
%     for l=(k-1)*l_time+1:k*l_time
%         Evec = -1*ones(l_space, 2*l_time);
%         for ts=1:2*l_time
%             Evec(:, ts) = V_wave_1(l_space*(ts-1)+1:l_space*ts,l);
%         end
%         subplot(2, ceil(set_size/2), it);
%         mesh(x_sm, y_sm, real(Evec'));
%         title(['EigVal ' num2str(D_zero(l,l))]);
%         xlabel('space');
%         ylabel('time');
%         it = it+1;
%     end
%     %filename = ['eigvec_wave_front_1_set_', num2str(k)];
%     %saveas(fig, filename);
% end


%%%%%%%%%%%%%% PLOT EIGENVECTORS w/ DOUBLE ENTRY %%%%%%%%%%%%%%

% how to plot eigenvectors, we have twice as many rows for time 



 
% [coeff_zero,score_zero] = pca(full(H_zero));
% [coeff_one,score_one] = pca(full(H_one));
% [coeff_wave_1,score_wave_1] = pca(full(H_wave_1));
% [coeff_wave_2,score_wave_2] = pca(full(H_wave_2));




