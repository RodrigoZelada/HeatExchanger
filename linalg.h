double det3x3(double A[3][3]){ 
   /*
    A = {A[0][0] A[0][1] A[0][2]
         A[1][0] A[1][1] A[1][2]
         A[2][0] A[2][1] A[2][2]
    }
    det = A[0][0]*|A[1][1] A[1][2]  - A[0][1]*|A[1][0] A[1][2]  +  A[0][2]*|A[1][0] A[1][1]
                   A[2][1] A[2][2]|            A[2][0] A[2][2] |            A[2][0]  A[2][1]   |  
         = A[0][0]*(A[1][1]*A[2][2] 
          - A[2][1]*A[1][2]) - A[0][1]*(A[1][0]*A[2][2] - A[2][0]*A[1][2])
          + A[0][2]*(A[1][0]*A[2][1] - A[2][0]*A[1][2])
    */
    double det = A[0][0]*(A[1][1]*A[2][2] - A[2][1]*A[1][2]) 
               -A[0][1]*(A[1][0]*A[2][2] - A[2][0]*A[1][2])
               +A[0][2]*(A[1][0]*A[2][1] - A[2][0]*A[1][1]);
    return det;
}

double det4x4(double A[4][4]){ 
   /*
    A = {A[0][0] A[0][1] A[0][2] A[0][3]
         A[1][0] A[1][1] A[1][2] A[1][3]
         A[2][0] A[2][1] A[2][2] A[2][3]
         A[3][0] A[3][1] A[3][2] A[3][3]
    }
    double A0[3][3] = {{A[1][1] A[1][2] A[1][3]},
          {A[2][1] A[2][2] A[2][3]},
          {A[3][1] A[3][2] A[3][3]}};
    double A1[3][3] = {{A[1][0] A[1][2] A[1][3]},
          {A[2][0] A[2][2] A[2][3]},
          {A[3][0] A[3][2] A[3][3]}};
    double A2[3][#] = {{A[1][0] A[1][1] A[1][3]},
            A[2][0] A[2][1] A[2][3]},
            A[3][0] A[3][1] A[3][3]}};
    double A3[3][3] = {{A[1][0] A[1][1] A[1][2] }.
            A[2][0] A[2][1] A[2][2] },
            A[3][0] A[3][1] A[3][2]}};
    double det = A[0][0]*det3x3(A0) - A[0][1]*det3x3(A1)+ A[0][2]*det3x3(A2) - A[0][3]*det3x3(A3);
    */
    double A0[3][3] = {{A[1][1], A[1][2], A[1][3]},
                    {A[2][1], A[2][2], A[2][3]},
                    {A[3][1], A[3][2], A[3][3]}};
    double A1[3][3] = {{A[1][0], A[1][2], A[1][3]},
                    {A[2][0], A[2][2], A[2][3]},
                    {A[3][0], A[3][2], A[3][3]}};
    double A2[3][3] = {{A[1][0], A[1][1], A[1][3]},
                    {A[2][0], A[2][1], A[2][3]},
                    {A[3][0], A[3][1], A[3][3]}};
    double A3[3][3] = {{A[1][0], A[1][1], A[1][2]},
                    {A[2][0], A[2][1], A[2][2] },
                    {A[3][0], A[3][1], A[3][2]}};
    double det = A[0][0]*det3x3(A0) - A[0][1]*det3x3(A1)+ A[0][2]*det3x3(A2) - A[0][3]*det3x3(A3);
    return det;
}

vector<double> cross(vector<double> AB, vector<double> AC){
      //b x c = {b1*c2 - b2*c1, b2*c0 - b0*c2 , b0*c1 - b1*c0}
	vector<double> BC = {AB[1]*AC[2] - AB[2]*AC[1], AB[2]*AC[0] - AB[0]*AC[2], AB[0]*AC[1] - AB[1]*AC[0]};
      return BC;
};


double dotArray(double a[],double b[], int N){
      double dot = 0.;
      for (int i=0;i<N;i++){
            dot += a[i]*b[i];
      }
      return dot;
}

double dot(vector<double> a, vector<double> b){
      int N = a.size();
      double dot = 0.;
      for (int i=0;i<N;i++){
            dot += a[i]*b[i];
      }
      return dot;
}

double norm(vector<double> a){
      int N = a.size();
      double dot = 0.;
      for (int i=0;i<N;i++){
            dot += a[i]*a[i];
      }
      return sqrt(dot);
}

