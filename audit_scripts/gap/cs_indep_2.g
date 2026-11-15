F := FreeGroup("u","v","w","x","e1","e2","e3","e4","e5","e6","e7","e8","e9","e10","e11","e12","e13","e14","e15","e16","e17","e18","e19","e20","e21","e22","e23","e24");;
u:=F.1;; v:=F.2;; w:=F.3;; x:=F.4;;
e1:=F.5;;
e2:=F.6;;
e3:=F.7;;
e4:=F.8;;
e5:=F.9;;
e6:=F.10;;
e7:=F.11;;
e8:=F.12;;
e9:=F.13;;
e10:=F.14;;
e11:=F.15;;
e12:=F.16;;
e13:=F.17;;
e14:=F.18;;
e15:=F.19;;
e16:=F.20;;
e17:=F.21;;
e18:=F.22;;
e19:=F.23;;
e20:=F.24;;
e21:=F.25;;
e22:=F.26;;
e23:=F.27;;
e24:=F.28;;
R := [
 e2*e4,
 e3*e2*e1*(e2)^-1,
 e6*e7,
 e4^-1*e6*e5*(e6)^-1,
 e8*e1^-1,
 e9*e8*e6^-1*(e8)^-1,
 e9^-1*e10,
 e10*e11,
 e11*e2^-1,
 e13*e3^-1,
 e14*e13*e12*(e13)^-1,
 e14^-1*e15,
 e15*e16,
 e16*e7^-1,
 e17*e8^-1,
 e18*e17*e5^-1*(e17)^-1,
 e19*e18^-1,
 e10^-1*e19*e17^-1*(e19)^-1,
 e20*e12^-1,
 e21*e20*e11^-1*(e20)^-1,
 e21^-1*e15^-1,
 e22*e21^-1*e13^-1*(e21^-1)^-1,
 e24*e16^-1,
 e19^-1*e24*e23*(e24)^-1,
 e23^-1*e20^-1,
 e22^-1*e23^-1*e24^-1*(e23^-1)^-1
];;
m := function(i) local es; es := [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1]; return GeneratorsOfGroup(F)[i+4]^es[i]; end;;
Add(R, x^-1*m(1));;
# try each other meridian as the second generator:  u = m(j) * x^-1
tryj := function(j)
  local RR, G, iso, M, FM, gm, n, ordl, ord, rws, red, gg, one, res, i, allred;
  RR := Concatenation(R, [ u^-1*m(j)*x^-1, v^-1*x*u*x^-1, w^-1*x*v*x^-1 ]);
  G := F/RR;
  iso := IsomorphismFpMonoid(G); M := Range(iso); FM := FreeMonoidOfFpMonoid(M);
  gm := GeneratorsOfMonoid(FM); n := Length(gm);
  # generators in increasing order: u,U,v,V,w,W, x,X, then all edge generators (largest: to be eliminated)
  ord := BasicWreathProductOrdering(FamilyObj(gm[1]), gm);
  rws := KnuthBendixRewritingSystem(M, ord);
  MakeConfluent(rws);
  red := function(g) return ReducedForm(rws, UnderlyingElement(Image(iso, g))); end;
  gg := GeneratorsOfGroup(G); one := red(One(G));
  Print("j=", j, ": confluent with ", Length(Rules(rws)), " rules\n");
  Print("  [u,v]=1: ", red(Comm(gg[1],gg[2]))=one, "  [u,w]=1: ", red(Comm(gg[1],gg[3]))=one, "  [v,w]=1: ", red(Comm(gg[2],gg[3]))=one, "\n");
  Print("  x w x^-1 = ", red(gg[4]*gg[3]*gg[4]^-1), "\n");
  for i in [5..Length(gg)] do Print("  e", i-4, " = ", red(gg[i]), "\n"); od;
end;;
tryj(2);
