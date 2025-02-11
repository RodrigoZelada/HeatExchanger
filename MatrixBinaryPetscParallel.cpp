//g++ MatrixBinaryPetsc.cpp -o MatrixBinaryPetsc -std=c++14 -larmadillo
//./MatrixBinaryPetsc output/case3Vol
//time ./MatrixBinaryPetsc output/case3Vol -ksp_view_mat ascii:A.txt -mat_view draw -draw_save
//time mpiexec -n 8 ./MatrixBinaryPetsc output/case3Vol -pc_type lu -pc_factor_mat_solver_type superlu_dist -ksp_monitor
//time mpiexec -n 8 ./MatrixBinaryPetsc output/case3Vol -ksp_type gmres -pc_type bjacobi -sub_pc_type ilu -ksp_monitor
//time mpiexec -n 8 ./MatrixBinaryPetsc output/case3Vol -ksp_type gmres -pc_type asm -sub_pc_type ilu -ksp_monitor

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <set>

#include <petsc.h>

using namespace std;

string output;

const int region1 = 3; //hot
const int region2 = 2; //cold
const int label_interface = 10;
const double k1=4e-3;
const double k2=1e-3;
const double lambda=1e-2;
const double rhocp1 = 1.; //instead of 1e-6 in diffusivity
const double rhocp2 = 1.; //instead of 1e-6 in diffusivity
const double epsilon = 1e-2;
const double ks = 1e-2;// 111.*1e-6; //thermal diffusivity D or alpha of cupper;

double T1(double x, double y, double z) { return 10.;} 
double T2(double x, double y, double z) { return 0.;}  

#include "linalg.h"
#include "mesh3D.h"
#include "trunc.h"
#include "medit.h"
#include "fespace.h"
#include "data.h"
#include "newIndicesWithoutBoundaryConditions.h"
#include "nitscheDiscontinuousVentcell3DParallel.h"
#include "nitscheAdjointDiscontinuousVentcell3DParallel.h"

int main(int argc, char ** argv) {
	output = argv[1];
	/*Duplicate degrees at interface*/
	mesh3D Th(output+"/Th.mesh");

	vector<int> vectorGlobalTriangles1(Th.nt);
	vector<int> vectorGlobalVertices1(Th.nv);
	vector<int> vectorGlobalTriangles2(Th.nt);
	vector<int> vectorGlobalVertices2(Th.nv);
	vector<int> localToGlobalVertices1;
	vector<int> localToGlobalVertices2;

	mesh3D Th1 = trunc(Th,region1,vectorGlobalTriangles1,vectorGlobalVertices1,localToGlobalVertices1);
	mesh3D Th2 = trunc(Th,region2,vectorGlobalTriangles2,vectorGlobalVertices2,localToGlobalVertices2);

	int lDBC1[1] = {3};
	int nDBC1 = 1;
	DirichletBC DBC1(Th1,lDBC1,nDBC1);
	int lDBC2[1] = {1};
	int nDBC2 = 1;
	DirichletBC DBC2(Th2,lDBC2,nDBC2);

	//only for Th2	
	/*int count=0; //all the points that are in Th2, that belongs to a DBC in Th1
	for (int i=0; i<Th2.nv; i++){
		int i_global = localToGlobalVertices2[i]; 
		if ( (DBC2.vertices[i][0] >=0 ) && (vectorGlobalVertices1[i_global] >= 0)){ //Th1 and Th2 share the vertex i
			int i1 = vectorGlobalVertices1[i_global];
			if (DBC1.vertices[i1][0] == -1){ //if is Dirichlet in Th1
				DBC2.vertices[i][0] = -1; DBC2.vertices[i][1] = DBC1.vertices[i1][1];
				count++;
			}
			else {
				DBC2.vertices[i][0] -= count;
			}
		}
		else{
			if (DBC2.vertices[i][0] >=0 ){ 
				DBC2.vertices[i][0] -= count;
			}
		}
	}
	DBC2.nvDBC += count; 
	DBC2.nvNotDBC -= count;*/

	//Bijection between triangles at the interface: t <-> t_tilde //
	//index triangle in Th1 || index triangle in Th2

	map<int, int> Gamma1;
	map<int, int> Gamma2;

	int t1,t2;
	for (int t_global=0; t_global<Th.nt; t_global++){
		if (Th.triangles[t_global][Th.d]==label_interface){
			t1 = vectorGlobalTriangles1[t_global];
			t2 = vectorGlobalTriangles2[t_global];
			Gamma1.insert({t1, t2});
			Gamma2.insert({t2, t1});
		}
	}

	//Bijection dirichlet - Th1 - Th2

	vector<double> kappa = readSolution(output+"/kappa.sol");
	vector<double>  kappa1 = truncP1(kappa,vectorGlobalVertices1);
	vector<double>  kappa2 = truncP1(kappa,vectorGlobalVertices2);
	vector<vector<double>>  u1 = readSolutionMat(output+"/uxHot.sol",output+"/uyHot.sol",output+"/uzHot.sol");
	vector<vector<double>>  u2 = readSolutionMat(output+"/uxCold.sol",output+"/uyCold.sol",output+"/uzCold.sol");

	int Nv1 = DBC1.nvNotDBC;
	int Nv2 = DBC2.nvNotDBC;
	int N = Nv1+Nv2; 

	Mat A;
    Vec x,b;
    KSP ksp;
	PetscInt Istart,Iend,chunk;
	PetscInitialize(&argc,&argv,NULL,"Solve laplacian\n");
    PetscMPIInt    rank,size; 
    int *recvcount, *displacements;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank); MPI_Comm_size(PETSC_COMM_WORLD,&size);

    VecCreate(PETSC_COMM_WORLD,&b);
    VecSetSizes(b,PETSC_DECIDE,N);
    VecSetFromOptions(b);
    VecGetOwnershipRange(b,&Istart,&Iend);

    MatCreate(PETSC_COMM_WORLD,&A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N);
    MatSetFromOptions(A);
	MatSetUp(A); 
    MatGetOwnershipRange(A,&Istart,&Iend);
    chunk = Iend-Istart;
    recvcount = new int[size]; displacements = new int[size];
    for (int i=0; i<size; i++){recvcount[i]=chunk;displacements[i]=Istart;}
    MPI_Allgather(&chunk,1,MPI_INT,recvcount,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&Istart,1,MPI_INT,displacements,1,MPI_INT,MPI_COMM_WORLD);
        	
	nitscheDiscontinuousVentcell3D(A, b, Th1, kappa1, u1, Gamma1, Th2, DBC1, DBC2, Istart, Iend);
	nitscheDiscontinuousVentcell3D(A, b, Th2, kappa2, u2, Gamma2, Th1, DBC1, DBC2, Istart, Iend);

	VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

	KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetOperators(ksp,A,A);
    KSPSetFromOptions(ksp);

	VecDuplicate(b,&x);
    VecSet(x,0.0);
    KSPSolve(ksp,b,x);
    VecViewFromOptions(x,NULL,"-vec_view");

	//SaveSolution

	PetscScalar *Tdistributed;
        VecGetArray(x,&Tdistributed); 

	PetscScalar T[N];
	MPI_Gatherv(Tdistributed,chunk,MPI_DOUBLE,&T,recvcount,displacements,MPI_DOUBLE,0,MPI_COMM_WORLD);
        
        if (rank==0){
        vector<double> uh1 = P1Function(Th1, &T1); vector<double> uh2 = P1Function(Th2, &T2);
		vector<double> uh1Global = P1Function(Th, &T1); vector<double> uh2Global = P1Function(Th, &T2);
		for (int i=0;i<Th1.nv;i++){
			if (DBC1.vertices[i][0] >= 0){
				uh1[i] = T[DBC1.vertices[i][0]];
            	uh1Global[localToGlobalVertices1[i]] = uh1[i]; uh2Global[localToGlobalVertices1[i]] = uh1[i];
			}
		}
		for (int i=0;i<Th2.nv;i++){
			if (DBC2.vertices[i][0] >= 0){
				uh2[i] = T[Nv1+DBC2.vertices[i][0]]; 
            	uh1Global[localToGlobalVertices2[i]] = uh2[i]; uh2Global[localToGlobalVertices2[i]] = uh2[i];
			}
			/*if (DBC2.vertices[i][1] == 3){ //Dirichlet belonging to Th1
				uh2[i] = T1(Th2.vertices[i][0],Th2.vertices[i][1],Th2.vertices[i][2]);
            	uh1Global[localToGlobalVertices2[i]] = uh2[i]; uh2Global[localToGlobalVertices2[i]] = uh2[i];
			}*/
		}

        cout << "uh1,uh2 obtained" << endl;

	//Compute J = int2d(Thhot,4)(rhocp1*T*u*n) - int2d(Thcold,2)(rhocp2*T*u*n)
	//and since n = (0,0,1) on 4 Thhot, n=(0,1,0) on 2 Thcold, then
	//J = int2d(Thhot,4)(rhocp1*T*uz) - int2d(Thcold,2)(rhocp2*T*uy)
	//

		double J=0.;
		saveSolution(output+"/T1.sol",uh1);saveSolution(output+"/T2.sol",uh2);
        saveSolution(output+"/T1Global.sol",uh1Global);saveSolution(output+"/T2Global.sol",uh2Global);
		vector<vector<double>> coef1 = baseSpace(Th1);
		vector<vector<double>> coef2 = baseSpace(Th2);
		vector<double> coefii(3);
		vector<double> gradii(3), uii(3);
		for (int n=0; n<Th1.ntet; n++){
			double vol = Th1.volume[n];
			gradii[0]=0.;gradii[1]=0.;gradii[2]=0.;
			for (int p=0; p<=Th1.d; p++){
				int pp = Th1.tetrahedra[n][p];
				for (int k=0;k<Th1.d;k++){
					coefii[k] = coef1[n][4*p+k]; gradii[k] += uh1[pp]*coefii[k];
				}
			}
			for (int i=0; i<=Th1.d; i++){
				int ii = Th1.tetrahedra[n][i];
				for (int k=0;k<Th1.d;k++){uii[k] = u1[ii][k];}
				J += (vol/4.)*rhocp1*dot(gradii,uii);
			}
		}
		//cold part
		for (int n=0; n<Th2.ntet; n++){
			double vol = Th2.volume[n];
			gradii[0]=0.;gradii[1]=0.;gradii[2]=0.;
			for (int p=0; p<=Th2.d; p++){
				int pp = Th2.tetrahedra[n][p];
				for (int k=0;k<Th1.d;k++){
					coefii[k] = coef2[n][4*p+k]; gradii[k] += uh2[pp]*coefii[k];
				}
			}
			for (int i=0; i<=Th2.d; i++){
				int ii = Th2.tetrahedra[n][i];
				for (int k=0;k<Th2.d;k++){uii[k] = u2[ii][k];}
				J -= (vol/4.)*rhocp2*dot(gradii,uii);
			}
		}
		ofstream file(output+"/J.gp");
		file.precision(8);
		file << J << endl;
		file.close();
        cout << "J = " << J << endl;	
       }
    KSPDestroy(&ksp);
    MatDestroy(&A);
    VecDestroy(&b);
    VecDestroy(&x);
	//Adjoint
	Mat AA;
    Vec xA,bA;
    KSP kspA;

    VecCreate(PETSC_COMM_WORLD,&bA);
    VecSetSizes(bA,PETSC_DECIDE,N);
    VecSetFromOptions(bA);
    VecGetOwnershipRange(bA,&Istart,&Iend);

    MatCreate(PETSC_COMM_WORLD,&AA);
    MatSetSizes(AA,PETSC_DECIDE,PETSC_DECIDE,N,N);
    MatSetFromOptions(AA);
	MatSetUp(AA); //PetscCall(MatSeqAIJSetPreallocation(A, 30, NULL)); 
	MatGetOwnershipRange(AA,&Istart,&Iend);
	
	nitscheAdjointDiscontinuousVentcell3D(AA, bA, Th1, kappa1, u1, Gamma1, Th2, DBC1, DBC2, Istart, Iend);
	nitscheAdjointDiscontinuousVentcell3D(AA, bA, Th2, kappa2, u2, Gamma2, Th1, DBC1, DBC2, Istart, Iend);

	VecAssemblyBegin(bA);
    VecAssemblyEnd(bA);
    MatAssemblyBegin(AA,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(AA,MAT_FINAL_ASSEMBLY);

	KSPCreate(PETSC_COMM_WORLD,&kspA);
    KSPSetOperators(kspA,AA,AA);
    KSPSetFromOptions(kspA);

	VecDuplicate(bA,&xA);
    VecSet(xA,0.0);
	KSPSolve(kspA,bA,xA);
    VecViewFromOptions(xA,NULL,"-vec_view");

	PetscScalar *TAdistributed;
    VecGetArray(xA,&TAdistributed);

	PetscScalar TA[N];
	MPI_Gatherv(TAdistributed,chunk,MPI_DOUBLE,&TA,recvcount,displacements,MPI_DOUBLE,0,MPI_COMM_WORLD); 

	if (rank==0){
        vector<double> vh1 = P1Function(Th1, &T2); vector<double> vh2 = P1Function(Th2, &T2); //T2
        vector<double> vh1Global = P1Function(Th, &T2); vector<double> vh2Global = P1Function(Th, &T2);
		for (int i=0;i<Th1.nv;i++){
			if (DBC1.vertices[i][0] >= 0){
				vh1[i] = TA[DBC1.vertices[i][0]];
                vh1Global[localToGlobalVertices1[i]] = vh1[i]; vh2Global[localToGlobalVertices1[i]] = vh1[i];
			}
		}
		for (int i=0;i<Th2.nv;i++){
			if (DBC2.vertices[i][0] >= 0){
				vh2[i] = TA[Nv1+DBC2.vertices[i][0]]; 
                vh1Global[localToGlobalVertices2[i]] = vh2[i]; vh2Global[localToGlobalVertices2[i]] = vh2[i];
			}
		}
		saveSolution(output+"/R1.sol",vh1);saveSolution(output+"/R2.sol",vh2);
        saveSolution(output+"/R1Global.sol",vh1Global);saveSolution(output+"/R2Global.sol",vh2Global);
	}
	if (rank==0){cout << "vh1,vh2 obtained" << endl;}

    KSPDestroy(&kspA);
    MatDestroy(&AA);
    VecDestroy(&bA);
    VecDestroy(&xA);
    
    PetscFinalize();
	
}
