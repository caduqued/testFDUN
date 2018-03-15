%% Prueba funcion bidimensional
clear all;
close all;
clc;

prueba01;
prueba02;
prueba03;

figure;ax=gca;
ax.XScale='log';ax.YScale='log';
loglog(errorVec01(:,2),errorVec01(:,3),'ko-');
hold on; grid on;
loglog(errorVec02(:,2),errorVec02(:,3),'kp-');
loglog(errorVec03(:,2),errorVec03(:,3),'kd-.');

