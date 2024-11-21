function [ MSS, VAL ] = DTC( VAL, MSS, POH )

for i = 1:size( POH, 2 )

    % Find dimension(s) for trisection
    max_L = max( MSS.LL(:, POH(i)) );
    ls = find( MSS.LL(:, POH(i)) == max_L );
    lls = length( ls );

    if ( lls ~= 1 ) && ( ~VAL.condition )

        % This is the least split side with tiebreak
        [ ~, indexas ] = min( VAL.minas( ls ) );
        ls = ls(indexas);
        lls = 1;

    end

    % Split sides
    VAL.minas( ls ) = VAL.minas( ls ) + 1;

    % Initialize
    l = MSS.LL(:, POH(i)) * ones( 1, 2 * lls );
    c = MSS.CC(:, POH(i)) * ones( 1, 2 * lls );
    O = zeros(1, 2 * lls);

    % Calculate new points
    c(ls, 1:2:end) = c(ls, 1:2:end) - diag( repelem( max_L / 3, lls ) );
    c(ls, 2:2:end) = c(ls, 2:2:end) + diag( repelem( max_L / 3, lls ) );

    % Evaluate points at the objective function
    f = arrayfun( @(x) VAL.TPD( c(:, x) ), (1:2 * lls) );

    % Find partitioning order
    [ ~, order ] = sort( [min( f(1:2:end), f(2:2:end) )', ls], 1 );

    for j = 1:lls

        % Calculate hyper-rectangles side lengths
        l(ls(order(1:j, 1)), (j * 2 - 1)) = max_L / 3;
        l(ls(order(1:j, 1)), (j * 2)) = max_L / 3;

        % Set partitioning order
        O((j * 2 - 1):(j * 2)) = [(order(j, 1) * 2 - 1), (order(j, 1) * 2)];

    end

    % Calculate diameters
    s = arrayfun( @(x) VAL.DiamL( l(:, x) ), (1:2 * lls) );

    % Store \ update information
    index                   = [ POH(i), VAL.I + 1:VAL.I + length(f) ];
    MSS.FF(index(2:end))    = f(:, O);
    MSS.DD(index)           = [ s(end), s ];
    MSS.LL(:, index)        = [ l(:,end), l ];
    MSS.CC(:, index(2:end)) = c(:, O);
    VAL.I                   = VAL.I + length(f);

    % Break if function evaluations exceed the budget
    if ( VAL.Limit <= VAL.I )
        break;
    end

end

% Update globl values
New_samples = flip( [VAL.Fminindex, VAL.feval + 1:VAL.I] );
[ Fmin, min_index ] = min( MSS.FF(New_samples) );
if VAL.Fmin >= Fmin
    VAL.Fmin = Fmin;
    VAL.Fminindex = New_samples(min_index);
    VAL.Xmin = MSS.CC(:, VAL.Fminindex);
end
VAL.feval = VAL.I;


return