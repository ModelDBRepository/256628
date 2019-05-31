%trap_integ.m
%Arjun Balachandar 2016
%Numerical method to calculate double integral, using trapezoid method

function [volumes] = trap_integ(fun,x_domain,y_domain,FP_domain)
    volumes = zeros(max(FP_domain(:))+1);
    xn = length(x_domain);
    yn = length(y_domain);
    
    for i=1:xn-1
        x2 = x_domain(i+1);
        x1 = x_domain(i);
        for j=1:yn-1
            y2 = y_domain(j+1);
            y1 = y_domain(j);
            delta = ( ((x2-x1)/2)*((y2-y1)/2)*(fun(x2,y2) + fun(x1,y2) + fun(x2,y1) + fun(x1,y1)) );
            volumes(FP_domain(i,j)+1) = volumes(FP_domain(i,j)+1) + delta;
        end
    end
    
    volumes = volumes(:,1);
end