//define P1 given functions
vector<double> P1Function(mesh3D &Th, double (*f)(double, double, double)){
  vector<double> TExact(Th.nv);
  for (int i=0; i < Th.nv; i++){
    TExact[i] = f(Th.vertices[i][0],Th.vertices[i][1],Th.vertices[i][2]);
  }
  return TExact;
}
