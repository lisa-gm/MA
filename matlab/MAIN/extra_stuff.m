%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% EXTRA STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% convergence things from basic_monodomain



% to measure convergence, do list of indices
%i = [3, 4, 5, 6, 22, 25, 28]; j = [3, 50, 40, 30, 5, 20, 30];
i = 5; j = 8;
k = (2.*j-3)*Nx + i;
% x_star = [1, 1, 1, 1, 0, 0, 0]';
% 
u_k = zeros(length(i), max_iter);
u_k(:,1) = u(k);

% in loop
u_k(:,iter+1) = u(k);


u_k = u_k(:,1:iter);
conv_it = zeros(length(i), iter-1);

for it=1:iter-1
    conv_it(:,it) = norm(u_k(:, it) - u_k(:, it+1), 1);
end

conv_it;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%


it = 1; it_row = 1;

% while it<iter-1
%     while it_row<length(i)
%         it+1
%         if(norm((u_k(it_row,it) - u_k(it_row,it+1)), 1) < 10^-4)
%             it
%             u_k = u_k(:,1:it+1);
%             it = iter;
%             it_row = length(i);
%         end
%         it = it+1; it_row = it_row+1;
%     end
% end

% u_M = zeros(length(i), size(u_k, 2)-1);
% convergence rate: quadratic conv. to x^* if 
% lim k -> \inf |x_k+1 - x^*|/|x_k - x^*|^2 < M, for some M>0

% for it=1:size(u_k,2)-1
%     u_M(:,it) = abs(u_k(:,it+1)- x_star)./abs(u_k(:,it) - x_star).^2;
% end

% figure
% plot(u_M');
% legend('1','2','3','4','5','6','7');
% title('convergence');
% 
% figure
% plot(u_k');
% title('point iteration u_k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% EIGENFUNCTION PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure;
% title('one zone');
% for l=1:l_time*l_space
%     Evec = -1*ones(l_space, l_time);
%     %some reshaping necessary, (0,0) is south-west corner of domain
%     for ts=1:l_time
%         Evec(:, ts) = V_one(l_space*(ts-1)+1:l_space*ts,l);
%     end
%     hold on;
%     mesh(x_sm, y_sm, Evec');
%     xlabel('space');
%     ylabel('time');
% end
% 
% figure;
% title('wavefront area 1');
% for l=1:16
%     %some reshaping necessary, (0,0) is south-west corner of domain
%     Evec = reshape(fliplr(V_wave_1(:,l)'), [l_time, l_space]);
%     hold on;
%     mesh(x_sm, y_sm, Evec);
%     xlabel('space');
%     ylabel('time');
% end
% 
% figure;
% title('wavefront area 2');
% for l=1:16
%     %some reshaping necessary, (0,0) is south-west corner of domain
%     Evec = reshape(fliplr(V_wave_2(:,l)'), [l_time, l_space]);
%     hold on;
%     mesh(x_sm, y_sm, Evec);
%     xlabel('space');
%     ylabel('time');
% end
% 
% figure;
% title('zero zone');
% for l=1:16
%     %some reshaping necessary, (0,0) is south-west corner of domain
%     Evec = reshape(fliplr(V_zero(:,l)'), [l_time, l_space]);
%     hold on;
%     mesh(x_sm, y_sm, Evec);
%     xlabel('space');
%     ylabel('time');
% end

% only every other point .... 


% we have no_of_sets many sets and set_size many points in there
% plot them individually 

%color_list = [ winter, hsv, summer, autumn, spring ];

% figure;
% suptitle('wave front 1');    
% for l=1:l_space*l_time
%     Evec = -1*ones(l_space, l_time);
%     for ts=1:l_time
%         Evec(:, ts) = V_wave_1(l_space*(ts-1)+1:l_space*ts,l);
%     end
%     subplot(set_size, no_of_sets, l);
%     mesh(x_sm, y_sm, Evec');
%     title(['EigVal ' num2str(D_wave_1(l,l))]);
%     xlabel('space');
%     ylabel('time');
% end
% 
% figure;
% suptitle('wave front 2');    
% for l=1:l_space*l_time
%     Evec = -1*ones(l_space, l_time);
%     for ts=1:l_time
%         Evec(:, ts) = V_wave_2(l_space*(ts-1)+1:l_space*ts,l);
%     end
%     subplot(set_size, no_of_sets, l);
%     mesh(x_sm, y_sm, Evec');
%     title(['EigVal ' num2str(D_wave_2(l,l))]);
%     xlabel('space');
%     ylabel('time');
% end
% 
% figure;
% suptitle('one zone');    
% for l=1:l_space*l_time
%     Evec = -1*ones(l_space, l_time);
%     for ts=1:l_time
%         Evec(:, ts) = V_one(l_space*(ts-1)+1:l_space*ts,l);
%     end
%     subplot(set_size, no_of_sets, l);
%     mesh(x_sm, y_sm, Evec');
%     title(['EigVal ' num2str(D_one(l,l))]);
%     xlabel('space');
%     ylabel('time');
% end