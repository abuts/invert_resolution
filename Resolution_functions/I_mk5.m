function int = I_mk5(m_v,n_t,v_avrg)

AbsErr = 0;
RelTol = 1.e-8;
Nu_v = v_avrg*v_avrg-0.25;

if m_v~=0
    if n_t==0
        int = 0; % not obvious, but the row below is always = 0
    else  
        num1 = 2i*pi*m_v;
        num2 = 2i*pi*n_t*Nu_v;
        funD = @(u)(exp(num1*u-num2./(v_avrg+u)));
        int = integral(funD,-0.5,0.5,'RelTol',0,'AbsTol',AbsErr);
    end
else
    if n_t == 0
        int = 1;
    else
        num = (-2i*pi*n_t*Nu_v);
        funI = @(u)(exp(num./(u+v_avrg)));
        
        int = integral(funI,-0.5,0.5,'RelTol',RelTol,'AbsTol',AbsErr);        
    end
end

