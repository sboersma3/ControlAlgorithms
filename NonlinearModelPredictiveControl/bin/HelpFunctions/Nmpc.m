function [uopt,output] = Nmpc(r,x0,d,p,ops)

u0 = ops.u0'; 
u0 = u0(:);

[cost,g,lbu,ubu,lbg,ubg] = costfunction_nonlinearconstraints(x0,d,r,p,ops);

prob    = struct('f', cost, 'x', ops.Us, 'g', g);
solver  = casadi.nlpsol('solver', 'ipopt', prob,ops.opts);

% Solve the NLP
if isfield(ops,'lam_x0')
    output  = solver('x0', u0, 'lbx', lbu, 'ubx', ubu, ...
        'lbg', lbg, 'ubg', ubg,'lam_x0',ops.lam_x0,'lam_g0',ops.lam_g0);
else
    output  = solver('x0', u0, 'lbx', lbu, 'ubx', ubu, ...
        'lbg', lbg, 'ubg', ubg);
end

Uopt    = full(output.x);
uopt    = reshape(Uopt,ops.k1,ops.nu)';

    function [cost,g,lbu,ubu,lbg,ubg] = costfunction_nonlinearconstraints(x0,d,r,p,ops)
                
        temp  = ops.F('x0',x0,'Us',ops.Us,'ds',d,'rs',r);
        
        cost  = temp.J;
        c     = temp.c;
        ceq   = temp.ceq;
        
        lbu   = ops.lbu;
        ubu   = ops.ubu;
        
        g     = [c ceq ops.dUs'];
        lbg   = ops.lbg;
        ubg   = ops.ubg;
               
    end

end


