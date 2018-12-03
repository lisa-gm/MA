clear all;
close all;
beep off;

addpath('../src');
addpath('../src/quad_nodes');
addpath('../exe');

Nvec   = 31;
extend = 0;

u_fh = @(x) sin(pi*x);
%u_fh = @(x) exp(-(x-0.5).^2/0.1^2);

diff = zeros(4,Nvec);

for extend=0:4
    
    Nx = Nvec;
    xfine = linspace(0,1,Nx);
    xcoarse = xfine(1:2:Nx);
    
    [R, P] = get_pr_spatial(xfine, xcoarse, extend);
    
    ufine   = u_fh(xfine).';
    ucoarse = u_fh(xcoarse).';
    
    
    err(1,extend+1) = norm(ucoarse - R*ufine, inf);
    err(2,extend+1) = norm(ufine - P*ucoarse, inf);
    
    diff(extend+1,:) = fft( ufine - P*ucoarse );
end

fig = figure(1);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 16, 16]);
set(fig, 'Color','white');

semilogy(0:Nvec-1, abs(diff)); hold on;
legend('p=1','p=3','p=5','p=7','p=9');
ylim([1e-10 1e4]);
xlim([0 Nvec-1]);
xlabel('k');
ylabel('abs(x_k)');


% convrate = log10(err(2,2:end)./err(2,1:end-1))./log10(Nvec(1:end-1)./Nvec(2:end));
% 
% fprintf('Expected convergence rate: %2i \n', 2*extend+1);
% fprintf('\n');
% for nn=1:length(Nvec)-1
%     fprintf('Measured rate:             %5.2f \n', convrate(nn));
% end