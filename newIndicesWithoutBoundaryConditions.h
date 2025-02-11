class DirichletBC{ //for each domain Omega_i
public:
    int nvNotDBC;
    int nvDBC;
    int **vertices;

    DirichletBC(){};

    DirichletBC(mesh3D &Th,int* lDBC, int nDBC){
        int **verticesDirichletBC; //put [index,label]
        int *verticesNotDirichletBC;
        int label;
        //by triangles
        nvDBC = 0;
        set<int> I_D;
        for (int t = 0; t < Th.nt; t++){ //without duplicates
            for (int l = 0; l<nDBC; l++){
                label = lDBC[l];
                if (Th.triangles[t][Th.d] == label){
                    for (int i=0; i<Th.d; i++){
                        I_D.insert(Th.triangles[t][i]);
                    }
                }
            }
        }
        vector<int> vc(I_D.begin(), I_D.end());
        nvDBC=vc.size();
        verticesDirichletBC = new int* [nvDBC];
        for (int i=0;i<nvDBC;i++){
            verticesDirichletBC[i] = new int[2];
        }
        for (int i=0;i<nvDBC;i++){
            verticesDirichletBC[i][0] = vc[i];
            /*fill label*/
            for (int l = 0; l<nDBC; l++){
                label = lDBC[l];
                int t=0;
                while (t<Th.nt){
                    if ((Th.triangles[t][Th.d] == label) && ((Th.triangles[t][0] == vc[i]) || (Th.triangles[t][1] == vc[i]) || (Th.triangles[t][2] == vc[i])) ){
                        verticesDirichletBC[i][1] = label;  
                        t = Th.nt; 
                    }
                    t++;
                }
            }

        }
        nvNotDBC = Th.nv - nvDBC;
        verticesNotDirichletBC = new int[nvNotDBC];

        int count=0;
        for (int j=0; j<verticesDirichletBC[0][0]; j++){
            verticesNotDirichletBC[count] = j;
            count++;
        }
        for (int i=0; i<nvDBC-1; i++){
            int i1=verticesDirichletBC[i][0]; int i2=verticesDirichletBC[i+1][0];
            for (int j=i1+1; j<i2; j++){
                verticesNotDirichletBC[count] = j;
                count++;
            }
        }
        for (int j=verticesDirichletBC[nvDBC-1][0]+1; j<Th.nv; j++){
            verticesNotDirichletBC[count] = j;
            count++;
        }

        vertices = new int* [Th.nv];
        for (int i=0;i<Th.nv;i++){
            vertices[i] = new int[2];
        }
        for (int i=0; i<nvNotDBC; i++){
            vertices[verticesNotDirichletBC[i]][0] = i; vertices[verticesNotDirichletBC[i]][1] = 0;
        }
        for (int i=0; i<nvDBC; i++){
            vertices[verticesDirichletBC[i][0]][0] = -1; vertices[verticesDirichletBC[i][0]][1] = verticesDirichletBC[i][1];
        }

        for (int i=0;i<nvDBC;i++){
            delete[] verticesDirichletBC[i];
        }
		delete[] verticesNotDirichletBC;
    };

    ~DirichletBC() {
        for (int i=0;i<nvNotDBC+nvDBC;i++){
            delete[] vertices[i];
        }
    }
};
