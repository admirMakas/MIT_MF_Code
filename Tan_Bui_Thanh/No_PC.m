function [P,Ri] = No_PC(zin)
% Default no preconditioner
P = ones(length(zin),1); Ri = P;