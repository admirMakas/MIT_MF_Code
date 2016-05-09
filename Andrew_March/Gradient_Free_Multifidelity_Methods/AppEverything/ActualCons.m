% if delta is negative TR constraint is not returned
%
function [c,ceq]=ActualCons(x,constr,constrNoAppr,constr_param,fid)
evalin('base','count=count+1;')

ceq=[];
c=constr(x,fid(2),constr_param{:});
otherConstrs=isa(constrNoAppr, 'function_handle');
if(otherConstrs)
    [c2,ceq0]=constrNoAppr(x,constr_param{:});
    c=[c;c2];
    ceq=[ceq,ceq0];
end