void savemesh(mesh3D &Th,string name){
  ofstream file(name);
  file << "MeshVersionFormatted 2\n";
  file << "\n";
  file << "Dimension 3\n";
  file << "\n";
  file << "Vertices\n";
  file << Th.nv << "\n";
  file.precision(15);
  for (int i=0; i<Th.nv; i++){
    file << Th.vertices[i][0] << " " << Th.vertices[i][1] << " " << Th.vertices[i][2] << " " << Th.vertices[i][3] << "\n";
  }
  file << "\n";
  file << "Tetrahedra\n";
  file << Th.ntet << "\n";
  for (int i=0; i<Th.ntet; i++){
    file << Th.tetrahedra[i][0]+1 << " " << Th.tetrahedra[i][1]+1 << " " << Th.tetrahedra[i][2]+1 << " " << Th.tetrahedra[i][3]+1 << " " << Th.tetrahedra[i][4] << "\n";
  }
  file << "\n";
  file << "Triangles\n";
  file << Th.nt << "\n";
  for (int i=0; i<Th.nt; i++){
    file << Th.triangles[i][0]+1 << " " << Th.triangles[i][1]+1 << " " << Th.triangles[i][2]+1 << " " << Th.triangles[i][3] << "\n";
  }
  file << "\n" ;
  file << "End";
  file.close();
}

void saveSolution(string name, vector<double> &uh){
  int N = uh.size();
  ofstream solution(name);
  solution << "MeshVersionFormatted 1\n";
  solution << "\n";
  solution << "Dimension 3\n";
  solution << "\n";
  solution << "SolAtVertices\n";
  solution << N << "\n";
  solution << "1 1 \n"; 
  solution.precision(15);
  for (int i=0; i<N; i++){
    solution << scientific << uh[i] << "\n";
  }
  solution << "\n" ;
  solution << "End";
  solution.close();
}

void saveSolutionGp(string name, vector<double> &uh){
  int N = uh.size();
  ofstream solution(name);
  solution << N << "\n";
  for (int i=0; i<N; i++){
    solution << uh[i] << "\n";
  }
  solution.close();
}

vector<double> readSolution(string name){
  ifstream file(name, ios::in);
  string line;
  int N;

  while (getline(file, line)) {
    if (line == "SolAtVertices") {
      getline(file, line);
      N = stoi(line);
      break;
    }
  }
  file.close();
  vector<double> uh(N);

  ifstream file2(name, ios::in);
  while (getline(file2, line)) {
    if (line == "1 1 ") {
      for (int i = 0; i < N; i++) {
        file2 >> uh[i];
      }
      break;
    }
  }
  file2.close();
  return uh;
}

vector<vector<double>> readSolutionMat(string name1, string name2, string name3){
  ifstream file(name1, ios::in);
  string line;
  int N;

  while (getline(file, line)) {
    if (line == "SolAtVertices") {
      getline(file, line);
      N = stoi(line);
      break;
    }
  }
  file.close();
  vector<vector<double>> uh(N,vector<double>(3));

  ifstream file1(name1, ios::in);
  while (getline(file1, line)) {
    if (line == "1 1 ") {
      for (int i = 0; i < N; i++) {
        file1 >> uh[i][0];
      }
      break;
    }
  }
  file1.close();

  ifstream file2(name2, ios::in);
  while (getline(file2, line)) {
    if (line == "1 1 ") {
      for (int i = 0; i < N; i++) {
        file2 >> uh[i][1];
      }
      break;
    }
  }
  file2.close();

  ifstream file3(name3, ios::in);
  while (getline(file3, line)) {
    if (line == "1 1 ") {
      for (int i = 0; i < N; i++) {
        file3 >> uh[i][2];
      }
      break;
    }
  }
  file3.close();
  return uh;
}

