clear all;
close all;
beep off;

addpath('../src');
addpath('../src/quad_nodes');
type = 'gauss-lobatto';
M = 4;
colls = collocation_sdc(0.0, 1.0, M, type);

filename = strcat('weights_', type, '_M', num2str(M), '.txt');
dlmwrite(filename, colls.Smat, 'delimiter','\t','precision', 18);
dlmwrite(filename, colls.nodes, '-append', 'delimiter','\t','precision', 18, 'roffset',1);
dlmwrite(filename, colls.weights, '-append', 'delimiter','\t','precision', 18, 'roffset',1);