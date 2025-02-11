vector<vector<double>> baseSpace(mesh3D &Th){
	vector<vector<double>> coef(Th.ntet,vector<double>((Th.d+1)*(Th.d+1)));
	double coefii[4];
	double A[4][4];
	double A0[4][4],A1[4][4],A2[4][4],A3[4][4];
	double vol=1.;

	for (int n = 0; n < Th.ntet; n++){
		for (int i = 0; i <= Th.d; i++){
			/*Aii = { {Th.verticesInTriangle[i][0][n],    Th.verticesInTriangle[i][1][n],1},
							{Th.verticesInTriangle[l(0)][0][n], Th.verticesInTriangle[l(0)][1][n],1},
							{Th.verticesInTriangle[l(1)][0][n], Th.verticesInTriangle[l(1)][1][n],1} };*/

			for (int j=0; j<= Th.d; j++){
				A[j][Th.d] = 1.;
				A0[j][Th.d] = 1.; A1[j][Th.d] = 1.; A2[j][Th.d] = 1.; A3[j][Th.d] = 1.;
				for (int k=0;k<Th.d;k++){
					A[j][k] = Th.vertices[Th.tetrahedra[n][(i+j) % (Th.d+1)]][k];
					A0[j][k] = Th.vertices[Th.tetrahedra[n][(i+j) % (Th.d+1)]][k]; A1[j][k] = Th.vertices[Th.tetrahedra[n][(i+j) % (Th.d+1)]][k]; A2[j][k] = Th.vertices[Th.tetrahedra[n][(i+j) % (Th.d+1)]][k]; A3[j][k] = Th.vertices[Th.tetrahedra[n][(i+j) % (Th.d+1)]][k];
				}
			}
			vol = det4x4(A); //6.*Th.volume[n];
			A0[0][0] = 1.;  //A0[1][0] = 0.; A0[2][0] = 0.;  A0[3][0] = 0.;  coefii[0] = det4x4(A0)/vol; 
			A1[0][1] = 1.;  //A1[1][1] = 0.; A1[2][1] = 0.;  A1[3][1] = 0.;  coefii[1] = det4x4(A1)/vol; 
			A2[0][2] = 1.;  //A2[1][2] = 0.; A2[2][2] = 0.;  A2[3][2] = 0.;  coefii[2] = det4x4(A2)/vol; 
			A3[0][3] = 1.;  //A3[1][3] = 0.; A3[2][3] = 0.;  A3[3][3] = 0.;  coefii[3] = det4x4(A3)/vol;
			for (int j=1; j <=Th.d; j++){
				A0[j][0] = 0.; A1[j][1] = 0.; A2[j][2] = 0.; A3[j][3] = 0.;
			}
			coefii[0] = det4x4(A0)/vol; coefii[1] = det4x4(A1)/vol; coefii[2] = det4x4(A2)/vol; coefii[3] = det4x4(A3)/vol; 
			for (int j=0; j<=Th.d; j++){
				coef[n][4*i+j] = coefii[j];
			} 
		}
	}
	return coef;
}
