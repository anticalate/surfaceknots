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
tryj(JJ);
