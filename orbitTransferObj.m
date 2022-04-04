function J = orbitTransferObj(z,~)
%maximize fuel
% m = z(5);
% J = -m;
%minimize time
tf = z(end);
J = tf;