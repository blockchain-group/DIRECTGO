function y = Michalewicz(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Michalewicz.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2376.htm
%
% Globally optimal solution:
%   if n = 1; f = -0.801303410098553; x(1) = [2.20290552094332];
%   if n = 2; f = -1.80130341009855;  x(2) = [1.57079632679490];
%   if n = 3; f = -2.76039467999456;  x(3) = [1.28499157179788];
%   if n = 4; f = -3.69885709846664;  x(4) = [1.92305846996689];
%   if n = 5; f = -4.687658179088148; x(5) = [1.72046977393433]; 
%   if n = 6; f = -5.687658179088148;  x(6) = [1.57079632679490];
%   if n = 7; f = -6.680885314444030;  x(7) = [1.45441397073036];
%   if n = 8; f = -7.663757350716241;  x(8) = [1.75608652112611];
%   if n = 9; f = -8.660151715641344;  x(9) = [1.65571741732520];
%   if n = 10; f = -9.660151715641344;  x(10) = [1.57079632679490];
%   if n = 11; f = -10.657482257192113;  x(11) = [1.49772880331924];
%   if n = 12; f = -11.649574998714796;  x(12) = [1.69661630056355];
%   if n = 13; f = -12.647817985597966;  x(13) = [1.63007608052931];
%   if n = 14; f = -13.647817985597966;  x(14) = [1.57079632679490];
%   if n = 15; f = -14.646400190319412;  x(15) = [1.51754611465925];
%   if n = 16; f = -15.641864818949971;  x(16) = [1.66606451159373];
%   if n = 17; f = -16.6408282327948;  x(17) = [1.61632864033842];
%   if n = 18; f = -17.6408282327948;  x(18) = [1.57079632679490];
%   if n = 19; f = -18.6399508750239;  x(19) = [1.52890700215990];
%   if n = 20; f = -19.6370135993494;  x(20) = [1.64745635802859];
%
% Variable bounds:
%   0 <= x(i) <= pi, i = 1...n
%   
% Problem Properties:
%   n  = any dimension;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 0;
    y.ng = 0;
    y.nh = 0;
    y.xl = @(i) 0;
    y.xu = @(i) pi;
    y.fmin = @(i) michalewicz_minimum(i);
    y.xmin = @(i) michalewicz_xmin(i);
    return
end
n = length(x);
m = 10;
s = 0;
for i = 1:n
    s = s + sin(x(i))*(sin(i*x(i)^2/pi))^(2*m);
end
y = -s;
end
%--------------------------------------------------------------------------
% Compute numerically the minimum of the Michalewicz function for a given
% dimension.
%--------------------------------------------------------------------------
function fmin = michalewicz_minimum(dim)
    fmin = 0;
    xmin = zeros(1,dim);
% compute the minimum for each dimension
    for d = 1:dim
% compute the location of the peak, which is very close to the minimum
        n = round(0.25*d -0.5);
        fraction = sqrt((2*n+1)/(2*d));
        
% if the fraction equals 0.5, the peak is located at the minimum.
        if fraction == 0.5
            fmin = fmin  + -1;
            xmin(d) = 0.5*pi;
            continue;
        end
        
% determine the search domain for ternary search
        if (fraction < 0.5)
            x0 = fraction*pi;
            x3 = 0.5*pi;
        else
            x0 = 0.5*pi;
            x3 = fraction*pi;
        end
        
% ternary search
        while(abs(x3-x0) > 1e-14)
            x1 = x0+(x3-x0)/3;
            f1 = -sin(x1).*(sin(d*x1.^2/pi).^20);
            
            x2 = x3-(x3-x0)/3;
            f2 = -sin(x2).*(sin(d*x2.^2/pi).^20);
            
            % update the search range
            if( f2 < f1 )
                x0 = x1;
            else
                x3 = x2;
            end
        end
        
        % update the values of the minimum
        xmin(d) = (x3+x0)/2;
        f = -sin(xmin(d)).*(sin(d*xmin(d).^2/pi).^20);
        fmin = fmin + f;
    end
end  
%--------------------------------------------------------------------------
% Compute numerically the minimum of the Michalewicz function for a given
% dimension.
%--------------------------------------------------------------------------
function xmin = michalewicz_xmin(dim)
    fmin = 0;
    xmin = zeros(1,dim);
% compute the minimum for each dimension
    for d = 1:dim
% compute the location of the peak, which is very close to the minimum
        n = round(0.25*d -0.5);
        fraction = sqrt((2*n+1)/(2*d));
        
% if the fraction equals 0.5, the peak is located at the minimum.
        if fraction == 0.5
            fmin = fmin  + -1;
            xmin(d) = 0.5*pi;
            continue;
        end
        
% determine the search domain for ternary search
        if (fraction < 0.5)
            x0 = fraction*pi;
            x3 = 0.5*pi;
        else
            x0 = 0.5*pi;
            x3 = fraction*pi;
        end
        
% ternary search
        while(abs(x3-x0) > 1e-14)
            x1 = x0+(x3-x0)/3;
            f1 = -sin(x1).*(sin(d*x1.^2/pi).^20);
            
            x2 = x3-(x3-x0)/3;
            f2 = -sin(x2).*(sin(d*x2.^2/pi).^20);
            
            % update the search range
            if( f2 < f1 )
                x0 = x1;
            else
                x3 = x2;
            end
        end
        
        % update the values of the minimum
        xmin(d) = (x3+x0)/2;
        f = -sin(xmin(d)).*(sin(d*xmin(d).^2/pi).^20);
        fmin = fmin + f;
    end
    xmin = xmin(dim);
end  