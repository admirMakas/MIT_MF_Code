function px = proj(x,kku,kkl,option)
% PROJ
% projection onto the strictly feasible set
% Copyright (c) 2007 by Tan Bui
% All Rights Reserved. Massachusetts Institute of Technology

meps = 0;
px=min(kku,x);
px=max(kkl,px);