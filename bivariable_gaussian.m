%Arjun Balachandar 2016
%Bivariate Normal Distribution
%Function defining a bivariate normal distribution (BND) with set describing
%parameters, used by fit_bivariate and fit_bivariate_sigma

function [fun] = bivariable_gaussian(p,x0,y0,sigma_x,sigma_y)
    A = 1/(2*pi*sigma_x*sigma_y*(1-p^2)^0.5);

    a = 1/(sigma_x^2);
    b = p/(sigma_x*sigma_y);
    c = 1/(sigma_y^2);

    %BND equation:
    fun = @(x,y) A*exp( - (a*(x-x0).^2 - 2*b*(x-x0).*(y-y0) + c*(y-y0).^2)/(2*(1-p^2)) );
end