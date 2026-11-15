# G = <x1,x2 | 4 relators>.  New generators: x = x1, u = x2*x1^-1, v = x u x^-1, w = x v x^-1.
# Monoid presentation on x,X,u,U,v,V,w,W with a wreath-product ordering in which x, X are largest.
F := FreeGroup("x","u","v","w");; x:=F.1;; u:=F.2;; v:=F.3;; w:=F.4;;
x1 := x;; x2 := u*x;;
R := [
 x2^-2*x1*x2*x1*x2^-1*x1^-2*x2*x1*x2*x1^-1,
 x2^-1*x1^-1*x2^-1*x1^-1*x2*x1*x2*x1^-1*x2^-1*x1^3,
 x2^-2*x1^-1*x2*x1*x2*x1^-2*x2^-1*x1*x2*x1,
 x2^-1*x1^-1*x2*x1^-1*x2^-1*x1^-1*x2*x1*x2^-1*x1^3,
 v^-1*x*u*x^-1, w^-1*x*v*x^-1 ];;
G := F/R;;
Print("abelian invariants ", AbelianInvariants(G), "\n");
iso := IsomorphismFpMonoid(G);; M := Range(iso);;
gensM := GeneratorsOfMonoid(M);; Print(gensM, "\n");
FM := FreeMonoidOfFpMonoid(M);;
gm := GeneratorsOfMonoid(FM);;
ord := BasicWreathProductOrdering(FamilyObj(gm[1]), [gm[3],gm[4],gm[5],gm[6],gm[7],gm[8],gm[1],gm[2]]);;
rws := KnuthBendixRewritingSystem(M, ord);;
Print("start KB ", Length(Rules(rws)), " rules\n");
MakeConfluent(rws);
Print("confluent, rules: ", Length(Rules(rws)), "\n");
red := function(g) return ReducedForm(rws, UnderlyingElement(Image(iso, g))); end;
gg := GeneratorsOfGroup(G);; gx:=gg[1];; gu:=gg[2];; gv:=gg[3];; gw:=gg[4];;
one := red(One(G));
Print("[u,v] trivial: ", red(Comm(gu,gv)) = one, "\n");
Print("[u,w] trivial: ", red(Comm(gu,gw)) = one, "\n");
Print("[v,w] trivial: ", red(Comm(gv,gw)) = one, "\n");
z := gx*gw*gx^-1;;
Print("x w x^-1 reduces to ", red(z), "\n");
Print("x w x^-1 = w v^-2 u : ", red(z) = red(gw*gv^-2*gu), "   or = w^2 v^-1 u : ", red(z) = red(gw^2*gv^-1*gu), "\n");
for r in Rules(rws) do Print(r, "\n"); od;
