function y = Goldstein_and_Pricec(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Goldstein_and_Price.m
%
% Original source: 
% - Na, J., Lim, Y., & Han, C. (2017). A modified DIRECT algorithm 
%   for hidden constraints in an LNG process optimization. Energy, 
%   126, 488–500. https://doi.org/10.1016/j.energy.2017.03.047  
%
% Globally optimal solution:
%   f* = 3.53893589646755
%   x* = (0.0486962986001838, -0.984635261794115)  
%
% Constraints (including variable bounds):
%   g(1): -((x(1)-1)^2)-((x(2)-1)^2)+0.9 <= 0;
%   g(2): -((x(1)+1)^2)-((x(2)+1)^2)+1.1 <= 0;
%         -2 <= x(1) <= 2;
%         -2 <= x(2) <= 2;
%   
% Problem Properties:
%   n  = 2;
%   #g = 2;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = (1+(x(1)+x(2)+1).^2.*(19-14.*x(1)+3.*x(1).^2-14.*x(2)+6.*x(1).*x(2)+3.*x(2).^2)).*(30+(2.*x(1)-3.*x(2)).^2.*(18-32.*x(1)+12.*x(1).^2+48.*x(2)-36.*x(1).*x(2)+27.*x(2).^2));
end