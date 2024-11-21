function Simpl = VertexTriangulation(n)

for i = 1:n + 1
    for j = 1:n
      if (i > j)
      	Simplex.P(i, j) = 1;
      else
      	Simplex.P(i, j) = 0;
      end
    end
end
k = 1;

Simpl(k).V = Simplex.P';

vt = 1:(n + 1);
finish = 1; % true
while (finish == 1)
  	eil = n;
  	check = 0; % false
  	while((check == 0) && (eil ~= 1))
        if (vt(eil) < n + 1)
            if(Simplex.P(eil, vt(eil)) == 0)
                check = 1; % true
            else
                vt(eil) = vt(eil) + 1;
            end
        else
            vt(eil) = 1;
            eil = eil - 1;
        end
    end
    if(eil == 1)
        finish = 0; % false
    end
    while((eil < n + 1) && (finish))
        for i = 1:n
            if(i == vt(eil))
                if(Simplex.P(eil - 1, vt(eil)) ~= 1)
                    Simplex.P(eil, vt(eil)) = 1;
                else
                    Simplex.P(eil, i) = Simplex.P(eil - 1, i);
                    if(vt(eil) < n)
                        vt(eil) = vt(eil) + 1;
                    else
                        vt(eil) = 1;
                    end
                end
            else
                Simplex.P(eil, i)=Simplex.P(eil - 1, i);
            end
        end
        eil = eil + 1;
    end
    if (k < factorial(n))
        k = k + 1;
        Simpl(k).V = Simplex.P';
    end
end

return