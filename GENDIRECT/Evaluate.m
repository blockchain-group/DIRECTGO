function y = Evaluate(Problem, Point)

global Item_History Item_Evals Item_Mmax Item_xMin Item_Iter%#ok<*GVMIS>
y = feval(Problem.f, Point);
if isnan(y) || isinf(y)
    y = 10^300;
end
Item_Evals = Item_Evals + 1;

if y < Item_History(end, 3) && Item_Evals <= Item_Mmax
    Item_xMin = Point;
    Item_History(end + 1, :) = [Item_Iter, Item_Evals, y, toc];
elseif Item_Evals == Item_Mmax
    Item_History(end + 1, :) = [Item_Iter, Item_Evals, Item_History(end, 3), toc];
end

end
