mesh3D trunc(mesh3D &Th, int region, vector<int> &vectorGlobalTriangles, vector<int> &vectorGlobalVertices, vector<int> &localToGlobalVertices){
    //I need maps global to local index
    
    int d=3;
    int nv=0;
    int nt=0;
    int ntet=0;
    double** vertices;
	int** triangles;
	int** tetrahedra;

    for (int i=0;i<Th.nv;i++){vectorGlobalVertices[i]=-1;}

    for (int n = 0; n<Th.ntet; n++){
        if (Th.tetrahedra[n][4] == region){
            for (int i=0; i<4; i++){
                vectorGlobalVertices[Th.tetrahedra[n][i]] = Th.tetrahedra[n][i];
            }
            ntet++;
        }
    }

    for (int i=0;i<Th.nv;i++){
        if (vectorGlobalVertices[i]>=0){nv++;}
    }

    localToGlobalVertices.reserve(nv);
    nv=0;
    for (int i=0;i<Th.nv;i++){
        if (vectorGlobalVertices[i]>=0){
            localToGlobalVertices[nv]=i;
            vectorGlobalVertices[i]=nv;
            nv++;
        }
    }

    tetrahedra = new int* [ntet];
    for (int i = 0; i < ntet; i++) {
        tetrahedra[i] = new int[d + 2];
    }
    
    vertices = new double* [nv];
    for (int i = 0; i < nv; i++) {
        vertices[i] = new double[d+1];
    }

    ntet=0;
    int v;
    for (int n=0; n<Th.ntet; n++){
        if (Th.tetrahedra[n][4] == region){
            tetrahedra[ntet][4] = Th.tetrahedra[n][4];
            for (int i=0; i<4; i++){
                v =  vectorGlobalVertices[Th.tetrahedra[n][i]];
                tetrahedra[ntet][i] = v;
                for (int j=0; j<4; j++){
                    vertices[v][j] = Th.vertices[Th.tetrahedra[n][i]][j];
                }
            }
            ntet++;
        }
    }

    for (int t=0; t<Th.nt; t++){
        if ( (vectorGlobalVertices[Th.triangles[t][0]] >= 0) && (vectorGlobalVertices[Th.triangles[t][1]] >= 0) && (vectorGlobalVertices[Th.triangles[t][2]] >= 0) ){
            vectorGlobalTriangles[t]=nt;
            nt++;
        }
    }

    triangles = new int* [nt];
    for (int i = 0; i < nt; i++) {
        triangles[i] = new int[d + 1];
    }

    nt=0;
    for (int t=0; t<Th.nt; t++){
        if ( (vectorGlobalVertices[Th.triangles[t][0]] >= 0) && (vectorGlobalVertices[Th.triangles[t][1]] >= 0) && (vectorGlobalVertices[Th.triangles[t][2]] >= 0) ){
            triangles[nt][3] = Th.triangles[t][3];
            for (int i=0; i<3; i++){
                triangles[nt][i] = vectorGlobalVertices[Th.triangles[t][i]];
            }
            nt++;
        }
    }
    mesh3D Th1(d,nv,nt,ntet,vertices,triangles,tetrahedra);
    return Th1;
}

vector<double> truncP1(vector<double> kappa, vector<int> &localToGlobalVertices){
    int N=localToGlobalVertices.size();
    vector<double> uh(N);
    for (int i=0;i<N;i++){
        uh[i] = kappa[localToGlobalVertices[i]];
    }
    return uh;
}