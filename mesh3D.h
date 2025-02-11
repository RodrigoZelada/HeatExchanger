class mesh3D {
public:
	int d; //dimension
	int nv; //number of vertices
	int nt; //number of triangles
	int ntet; //number of tetrahedra

	double** vertices;
	int** triangles;
	int** tetrahedra;
	double* volume;
    double* area;
	double** normal;
    double* triangleToTetra;

	mesh3D(){};

	mesh3D(int d1,int nv1,int nt1,int ntet1,double** vertices1,int** triangles1,int** tetrahedra1)
	: d(d1),nv(nv1),nt(nt1),ntet(ntet1),vertices(vertices1),triangles(triangles1),tetrahedra(tetrahedra1) {
		trianToTetra();
		computeVolume();
		computeArea();
	}

	mesh3D(string meshFile) {
		loadMesh(meshFile);
		triangleToTetra = new double [nt];
		volume = new double [ntet];
		area = new double [nt];
		normal = new double* [nt];
		for (int i = 0; i < nt; i++) {
			normal[i] = new double[d];
		}
		/*trianToTetra();
		computeVolume();
		computeArea();*/
	};

	~mesh3D() {
		for (int i = 0; i < nv; i++) {
			delete[] vertices[i];
		}
		for (int i = 0; i < nt; i++) {
			delete[] triangles[i];
		}
		for (int i = 0; i < ntet; i++) {
			delete[] tetrahedra[i];
		}
		
		delete[] triangleToTetra;
		delete[] volume;
		delete[] area;
		for (int i = 0; i < nt; i++) {
			delete[] normal[i];
		}
		
	}

	template <typename T>
	void readField(ifstream &file, T **Array, int n, int cols, T dif) {
		T tmp;
		string line;
		for (int i = 0; i < n; i++) {
			getline(file, line);
			istringstream iss(line);

			int j = 0;

			while (iss >> tmp) {
				if (j<cols-1){
					Array[i][j] = tmp-dif;
				}
				else {
					Array[i][j] = tmp;
				}

				j++;
			}
		}
	}

	void loadMesh(string meshFile) {

		ifstream file(meshFile, ios::in);
		string line;

		while (getline(file, line)) {
			if (line == "Dimension") {
				getline(file, line);
				d = stoi(line);
			}
			if (line == "Dimension 3") {
				d = 3;
			}

			if (line == "Vertices") {
				getline(file, line);
				nv = stoi(line);

				vertices = new double* [nv];
				for (int i = 0; i < nv; i++) {
					vertices[i] = new double[d+1];
				}

				readField<double>(file, vertices, nv, d+1, 0.);
			}
			if (line == "Triangles") {
				getline(file, line);
				nt = stoi(line);

				triangles = new int* [nt];
				for (int i = 0; i < nt; i++) {
					triangles[i] = new int[d + 1];
				}

				readField<int>(file, triangles, nt, d+1, 1);
			}
			if (line == "Tetrahedra") {
				getline(file, line);
				ntet = stoi(line);

				tetrahedra = new int* [ntet];
				for (int i = 0; i < ntet; i++) {
					tetrahedra[i] = new int[d + 2];
				}

				readField<int>(file, tetrahedra, ntet, d+2, 1);
			}
		}
	}

 void computeVolume(){
	volume = new double [ntet];

    /*Compute volume*: two options: 1. Using dot and cross product or 2. determinant*/
    //1. By using |a dot (b x c) |
    vector<double> A(d),B(d),C(d),D(d),AB(d),AC(d),AD(d),BC(d);
    double vol;

    //double Pe[4][4];
    //double det,det1,det2,det3,det4;
    for (int n = 0; n < ntet; n++){
        for (int i=0;i<d;i++){
            A[i] = vertices[tetrahedra[n][0]][i];
            B[i] = vertices[tetrahedra[n][1]][i];
            C[i] = vertices[tetrahedra[n][2]][i];
            D[i] = vertices[tetrahedra[n][3]][i];

            AB[i] = B[i] - A[i];
            AC[i] = C[i] - A[i];
            AD[i] = D[i] - A[i];
        }
        //b x c = {b1*c2 - b2*c1, b2*c0 - b0*c2 , b0*c1 - b1*c0}
		BC = cross(AB,AC);
        vol = dot(BC,AD); 
        volume[n] = (1./6.)*fabs(vol);

        /*
        for (int k=0;k<=Th.d;k++){
            Pe[0][k] = 1.;
        }
        for (int j=1;j<=Th.d;j++){
            for (int k=0;k<=Th.d;k++){
                Pe[j][k] = Th.vertices[Th.tetrahedra[n][k]][j];
            }
        }
        vol = det4x4(Pe);
        Th.volume[n] = (1./6.)*fabs(vol);*/
    }

  };

void computeArea(){
	area = new double [nt];
	normal = new double* [nt];
	for (int i = 0; i < nt; i++) {
		normal[i] = new double[d];
	}

    double normVal, dotVal;

	int tetra;
	int i0,i1,i2,j;
	vector<double> A(d),B(d),C(d),D(d),AB(d),AC(d),AD(d),BC(d),P(d),x(d),N(d);
	j=0;
    for (int n = 0; n < nt; n++){
		tetra = triangleToTetra[n];
		i0 = triangles[n][0]; i1 = triangles[n][1]; i2 = triangles[n][2];
		for (int p=0; p<=d;p++){
			if ( (tetrahedra[tetra][p] != i0) && (tetrahedra[tetra][p] != i1) && (tetrahedra[tetra][p] != i2) ){
				j = tetrahedra[tetra][p];
			}
		}

        for (int i=0;i<d;i++){
		    A[i] = vertices[i0][i];  
            B[i] = vertices[i1][i];  
            C[i] = vertices[i2][i];  
            P[i] = vertices[j][i];  

            AB[i] = B[i] - A[i];
            AC[i] = C[i] - A[i];
        }
        //a x b = {a1*b2 - a2*b1,a2*b0 - a0*b2 ,a0*b1 - a1*b0}
		N = cross(AB,AC);
		normVal = norm(N);
        for (int i=0;i<d;i++){
            N[i] = N[i]/normVal; 
        }
		area[n] = 0.5*normVal;
        for (int i=0;i<d;i++){
		    x[i] = A[i] - P[i];  
        }
        dotVal = dot(x,N);
		if (dotVal < 0){
			normal[n][0] = -N[0];normal[n][1] = -N[1];normal[n][2] = -N[2];
		}
		else{
	        normal[n][0] = N[0];normal[n][1] = N[1];normal[n][2] = N[2];
		}
	}
  };

void trianToTetra(){
    triangleToTetra = new double [nt];

    for (int t=0; t<nt; t++){
      int vt0 = triangles[t][0]; int vt1 = triangles[t][1]; int vt2 = triangles[t][2];
      int tet=0;
      while (tet<ntet){
        int v0 = tetrahedra[tet][0]; int v1 = tetrahedra[tet][1]; int v2 = tetrahedra[tet][2]; int v3 = tetrahedra[tet][3];

        if ( ( (v0==vt0) || (v1==vt0) || (v2==vt0) || (v3==vt0) ) && ((v0==vt1) || (v1==vt1) || (v2==vt1) || (v3==vt1)) && ((v0==vt2) || (v1==vt2) || (v2==vt2) || (v3==vt2)) ){
          triangleToTetra[t] = tet;
          tet = ntet;
        }
        tet++;
      }
    }
  }

};
