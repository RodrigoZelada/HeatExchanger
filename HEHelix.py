from pyfreefem import FreeFemRunner, readFFArray
from nullspace_optimizer import Optimizable, nlspace_solve
from pymedit import Mesh, P1Function, trunc, Mesh3D, cube, mmg3d, generate3DMesh, P1Function3D, trunc3DMesh
from pymedit import saveToVtk, advect, P1Vector3D, mshdist
import numpy as np
import shutil
import pandas as pd

output='output/caseHelix' #1
import os
os.makedirs(output,exist_ok=True)

config={'OUTPUT':output}
#os.system("g++ IsolantComplete3DExt.cpp -o IsolantComplete3DExt -std=c++14 -larmadillo -fopenmp")

N=64 #64 #number of cores
i = format(0,'04d')
meshes='Helix'

r=0.5
xmax=-1
ymax=2.25
zmax=5.5
xmin=-1
ymin=-2.25
zmin=-5.5

barrier = lambda x,y,z: min(max((x-xmin)**2 + (y-ymin)**2 - r**2, z-zmin), max((x-xmax)**2 + (y-ymax)**2 - r**2,zmax-z) )
cutIn = lambda x,y,z: min((y-ymin)**2 + (x-xmin)**2 - r**2, zmin-z)
cutOut = lambda x,y,z: min((y-ymax)**2 + (x-xmax)**2 - r**2, z-zmax)
cut = lambda x,y,z: max(cutIn(x,y,z),cutOut(x,y,z))

#dimensions: [-3.5,3.5] x [-4.4,4.4] x [-7,7]
Th0 = Mesh3D(meshes+"/Helix.mesh")

# Meshing param
hmin = 1e-1 #2.5e-2 # 0.5
hmax = 3e-1 #5e-2 #1e-1 Maximum edge size
hgrad = 1.3 # Gradation factor
hausd = 5e-2 #1e-2 # Approximation quality factor, leave as it is

# Local parameter to prescribe local mesh size on the interface boundary 10 (corresponding to the level set)
refermmg = f"""Parameters
6

1 Triangles 6e-2 8e-2 5e-2
2 Triangles 6e-2 8e-2 5e-2
3 Triangles 6e-2 8e-2 5e-2
4 Triangles 6e-2 8e-2 5e-2
5 Triangles 8e-2 1e-1 5e-2
10 Triangles 6e-2 1e-1 5e-2
"""

Thinit = mmg3d(Th0,extra_args='-nomove -noswap -nosurf -noinsert',nr=False)
Thinit.save(output+'/Th0.mesh')

dOmega1 = mshdist(Thinit, ncpu=N)
dOmega1.save(output+'/d1.sol')
dOmega1.save(output+'/d1.o.sol')


preamble="""
func int readSolFile(mesh & Th, string fileName, real[int] & phi){
    ifstream f(fileName);
    string dummy="";
    while(dummy(0:2)!="Sol"){
        f>>dummy;
    }
    int n;
    f >> n;
    if(n!=Th.nv){
        cout << "Error : the number of vertices in the file "+fileName+" does not correspond to the mesh in memory"<< endl;
        exit(1);
    }
    f >> dummy;
    f >> dummy;
    for(int i=0;i<Th.nv;i++){
        f>>phi[i];
    }
}

func int saveArray(string fileName, real[int] &value){
    ofstream file(fileName);
    file.precision(16);
    file << value;
}

func int readData(string fileName, real[int] &data){
    {
        ifstream f(fileName);
        f>>data;
    }
}


include "macros.edp"
load "medit"
load "msh3"

int hotfluid = 3;
int coldfluid = 2;

real rhoCold=1., rhoHot=1.;
real nuCold=2., nuHot=1.;
real kHot=2e-3, kCold=1e-3, ks= 1e-2;
real cpCold=1e0, cpHot=1e0;
real maxVolHot = 96.4952019791, maxVolCold = 358.0544629510797;
real epsilon = 1e-2; 
"""

mesh_code='''

mesh3 ThBox=readmesh3("$MESH");
mesh3 Th=trunc(ThBox, ((region == coldfluid) || (region == hotfluid)));
savemesh(Th,"$OUTPUT/Th.mesh");
mesh3 Thcold = trunc(Th, (region==coldfluid));
mesh3 Thhot = trunc(Th, (region==hotfluid));
savemesh(Thcold,"$OUTPUT/Thcold.mesh");
savemesh(Thhot,"$OUTPUT/Thhot.mesh");
real Vcold = int3d(Thcold,qforder=1)(1.); real Vhot = int3d(Thhot,qforder=1)(1.);
cout << "Vcold = " << Vcold << ", Vhot = " << Vhot << endl; 
{
    ofstream f("$OUTPUT/Vcold.gp");
    f.precision(16);
    f << Vcold - maxVolCold  << endl;
}
{
    ofstream f("$OUTPUT/Vhot.gp");
    f.precision(16);
    f << Vhot - maxVolHot  << endl;
}

'''

solve_kappa = """
load "PETSc"
macro dimension()3 //EOM
include "macro_ddm.idp"

mesh3 Th=readmesh3("$OUTPUT/Th.mesh");
mesh3 ThGlobal = Th;

int[int] n2o;
macro ThN2O() n2o // this tells buildDmesh to keep the local to global correspondence
buildDmesh(Th)

fespace Gh(Th, [P1,P1,P1]);
fespace GhGlobal(ThGlobal, [P1,P1,P1]);
fespace Ph(Th, P1);
fespace PhGlobal(ThGlobal, P1);

Gh [gx, gy, gz];
GhGlobal [gxGlobal, gyGlobal, gzGlobal], [gxsum, gysum, gzsum];

Ph kappa, dOmega;
PhGlobal kappaGlobal, kappasum, dOmegaGlobal, nx, ny, nz;
dOmegaGlobal[] = readsol("$OUTPUT/d1.sol");

int[int] subIdx  = restrict(Ph, PhGlobal, n2o);

dOmega[] = dOmegaGlobal[](subIdx);

varf probgsd([gx,gy,gz],[thetax, thetay,thetaz]) = int3d(Th,qforder=3)( tr([gx,gy,gz])*[thetax, thetay, thetaz]      )
                                                   +int3d(Th,qforder=2)( tr(grad(dOmega))*[thetax, thetay, thetaz]  );

Mat AG; //Global distributed matrix
{
macro def(i)[i, i#B, i#C] //EOM
macro init(i)[i, i, i] // EOM
createMat(Th, AG, [P1,P1,P1])
}

AG=probgsd(Gh,Gh,tgv=-1);
set(AG, sparams = "-ksp_type cg -ksp_rtol 1e-6 -pc_type hypre"); //"-pc_type lu"
real[int] bG=probgsd(0,Gh,tgv=-1);
gx[] = AG^-1 * bG;

Mat A;
createMat(Th,A,P1)
varf probkappa(kappa, psi) = int3d(Th,qforder=3)(9*Th.hmax*Th.hmax*tr(grad(kappa))*grad(psi) + kappa*psi)
                             +int3d(Th,qforder=2)(div(gx,gy,gz)*psi);


A = probkappa(Ph, Ph, tgv=-1);
set(A, sparams = "-pc_type lu " );
real[int] b = probkappa(0, Ph, tgv=-1);
kappa[] = A^-1 * b;

kappa[] .*= A.D;
kappaGlobal[](subIdx) = kappa[]; //local function from #0 interpolated by 0, on the global domain
mpiReduce(kappaGlobal[], kappasum[], processor(0, mpiCommWorld), mpiSUM); //problem because of global elements => need to use partition of unity

if (mpirank==0){
    //save the solution
    savesol("$OUTPUT/kappa.sol",ThGlobal,kappasum);
    savesol("$OUTPUT/signeddistance.sol",ThGlobal,dOmegaGlobal);
}
"""

def computeU(mesh):
    M = Mesh3D(mesh)
    FreeFemRunner([preamble,mesh_code],config,run_dir=output,run_file='mesh_code.edp',debug=1).execute({'MESH':mesh})
    os.system(f"ff-mpirun -np {N} HelixNSCold.edp -output {output} -v 0")
    os.system(f"ff-mpirun -np {N} HelixNSHot.edp -output {output} -v 0")
    FreeFemRunner([preamble, solve_kappa],config,run_dir=output,run_file='solve.edp',debug=1).execute({'MESH':mesh},ncpu=N)
    #os.system(f"time mpiexec -n 1 ./MatrixBinaryPetscHelix {output} -pc_type lu -pc_factor_mat_solver_type superlu_dist -ksp_monitor ")
    os.system(f"time mpiexec -n {N} ./MatrixBinaryPetscHelixParallel {output} -pc_type lu -pc_factor_mat_solver_type superlu_dist -ksp_monitor ")
    os.system(f"ff-mpirun -np {N} adjointNSCold.edp -output {output} -v 0")
    os.system(f"ff-mpirun -np {N} adjointNSHot.edp -output {output} -v 0")
    os.system(f"ff-mpirun -np {N} adjointNSDCold.edp -output {output} -v 0")
    os.system(f"ff-mpirun -np {N} adjointNSDHot.edp -output {output} -v 0")

    with open(output+'/J.gp','r') as f:
        J = float(f.readlines()[0])
    with open(output+'/Dcold.gp','r') as f:
        Dcold = float(f.readlines()[0])
    with open(output+'/Dhot.gp','r') as f:
        Dhot = float(f.readlines()[0])
    with open(output+'/Vhot.gp','r') as f:
        Vhot = float(f.readlines()[0])
    return (J,Dcold,Dhot,Vhot)

code_sensitivity='''
load "PETSc"
macro dimension()3 //EOM
include "macro_ddm.idp"

mesh3 ThBox=readmesh3("$MESH");
mesh3 Th=ThBox;
mesh3 Thcold = readmesh3("$OUTPUT/Thcold.mesh");
mesh3 Thhot = readmesh3("$OUTPUT/Thhot.mesh");

int[int] n2o;
macro ThBoxN2O() n2o // this tells buildDmesh to keep the local to global correspondence

buildDmesh(ThBox)

fespace Ph(Th,P1);
Ph d1,T1Global, T2Global, R1Global, R2Global;
d1[]=readsol("$OUTPUT/d1.sol");

fespace Vhcold(Thcold, [P2,P2,P2,P1]); //global
Vhcold [uxCold, uyCold, uzCold, pCold], [vxCold, vyCold, vzCold, qCold], [vxDCold, vyDCold, vzDCold, qDCold];
readData("$OUTPUT/uxGlobalCold.gp",uxCold[]);
readData("$OUTPUT/vxGlobalCold.gp",vxCold[]);
readData("$OUTPUT/vxDGlobalCold.gp",vxDCold[]);

fespace Vhhot(Thhot, [P2,P2,P2,P1]); //global
Vhhot [uxHot, uyHot, uzHot, pHot], [vxHot, vyHot, vzHot, qHot], [vxDHot, vyDHot, vzDHot, qDHot];
readData("$OUTPUT/uxGlobalHot.gp",uxHot[]);
readData("$OUTPUT/vxGlobalHot.gp",vxHot[]);
readData("$OUTPUT/vxDGlobalHot.gp",vxDHot[]);

//Compute the shape derivative

fespace Vh(ThBox,[P1,P1,P1]); // space for deformation field theta
Vh [thetaxJ,thetayJ,thetazJ], [thetaxDCold,thetayDCold,thetazDCold], [thetaxDHot,thetayDHot,thetazDHot],  [thetaxVHot,thetayVHot,thetazVHot];

Vh [nx,ny,nz];
fespace PhParallel(ThBox,P1);
PhParallel kappa, dOmega;
int[int] subIdxK  = restrict(PhParallel, Ph, n2o);

dOmega[] = d1[](subIdxK);

varf probgsd([gx,gy,gz],[thetax, thetay,thetaz]) = int3d(ThBox,qforder=3)( tr([gx,gy,gz])*[thetax, thetay, thetaz]      )
                                                   +int3d(ThBox,qforder=2)( tr(grad(dOmega))*[thetax, thetay, thetaz]  );

Mat AG; //Global distributed matrix
{
macro def(i)[i, i#B, i#C] //EOM
macro init(i)[i, i, i] // EOM
createMat(ThBox, AG, [P1,P1,P1])
}

AG=probgsd(Vh,Vh,tgv=-1);
set(AG, sparams = "-ksp_type cg -pc_type hypre"); //"-pc_type lu"
real[int] bG=probgsd(0,Vh,tgv=-1);
nx[] = AG^-1 * bG;

Mat AK;
createMat(ThBox,AK,P1)
varf probkappa(kappa, psi) = int3d(ThBox,qforder=3)(9*ThBox.hmax*ThBox.hmax*tr(grad(kappa))*grad(psi) + kappa*psi)
                             +int3d(ThBox,qforder=2)(div(nx,ny,nz)*psi);


AK = probkappa(PhParallel, PhParallel, tgv=-1);
set(AK, sparams = "-pc_type lu " );
real[int] bK = probkappa(0, PhParallel, tgv=-1);
kappa[] = AK^-1 * bK;

fespace Phhot(Thhot,P1);
fespace Phcold(Thcold,P1);

Phhot THot,RHot;
Phcold TCold,RCold;

THot[]=readsol("$OUTPUT/T1.sol"); TCold[]=readsol("$OUTPUT/T2.sol"); T1Global[]=readsol("$OUTPUT/T1Global.sol"); T2Global[]=readsol("$OUTPUT/T2Global.sol"); 
RHot[]=readsol("$OUTPUT/R1.sol"); RCold[]=readsol("$OUTPUT/R2.sol"); R1Global[]=readsol("$OUTPUT/R1Global.sol"); R2Global[]=readsol("$OUTPUT/R2Global.sol");

fespace Ph2(Th,P2);
Ph2 f;
f = 0.5*(T1Global-T2Global)*(R1Global+R2Global);

varf shapeDerivativeJ([thetax,thetay,thetaz],[thetaxp,thetayp,thetazp])
    = int2d(ThBox,10,optimize=0)( ((2.*nuHot*tr(EPS(vxHot,vyHot,vzHot))*EPS(uxHot,uyHot,uzHot) - 2.*nuCold*tr(EPS(vxCold,vyCold,vzCold))*EPS(uxCold,uyCold,uzCold))  
     + (kCold*tr(grad(TCold))*grad(RCold) - kHot*tr(grad(THot))*grad(RHot)) 
    + 2.*(kHot*(dx(THot)*nx + dy(THot)*ny + dz(THot)*nz )*(dx(RHot)*nx + dy(RHot)*ny + dz(RHot)*nz ) - kCold*(dx(TCold)*nx + dy(TCold)*ny + dz(TCold)*nz )*(dx(RCold)*nx + dy(RCold)*ny + dz(RCold)*nz ) ) 
    -  epsilon*ks*kappa*tr(gradtaumeantauN(THot,TCold,nx,ny,nz))*gradtaumeantauN(RHot,RCold,nx,ny,nz)  +  epsilon*ks*kappa*2.*tr(MatrixByVector(GRAD(nx, ny, nz),gradtaumeantauN(THot,TCold,nx,ny,nz)))*gradtaumeantauN(RHot,RCold,nx,ny,nz)  - (ks/epsilon)*kappa*jumptau(TCold,THot)*jumptau(RCold,RHot)
    - ks*kappa*kappa*jumptau(TCold,THot)*meantau(RCold,RHot) + ks*tr(GRAD(nx,ny,nz))*GRAD(nx,ny,nz)*jumptau(TCold,THot)*meantau(RCold,RHot)) 
    * (thetaxp*nx + thetayp*ny + thetazp*nz))
    -int2d(ThBox,10,optimize=0)(ks*tr(gradtauN(f,nx,ny,nz))*MatrixByVector(GRADT(thetaxp,thetayp,thetazp),[nx,ny,nz]) +  ks*tr(gradtauN(f,nx,ny,nz))*MatrixByVector(GRADT(nx,ny,nz),[thetaxp,thetayp,thetazp])) 
    +on(1,2,3,4,5,6,thetax=0,thetay=0,thetaz=0);

varf shapeDerivativeVCold([thetax,thetay,thetaz],[thetaxp,thetayp,thetazp])
    =int2d(ThBox,10,optimize=0)(-thetaxp*nx - thetayp*ny - thetazp*nz)
     +on(1,2,3,4,5,6,thetax=0,thetay=0,thetaz=0);

varf shapeDerivativeVHot([thetax,thetay,thetaz],[thetaxp,thetayp,thetazp])
    =int2d(ThBox,10,optimize=0)(thetaxp*nx+thetayp*ny+thetazp*nz)
     +on(1,2,3,4,5,6,thetax=0,thetay=0,thetaz=0);

varf shapeDerivativeDCold([thetax,thetay,thetaz],[thetaxp,thetayp,thetazp])
    =int2d(ThBox,10,optimize=0)( 
                        (-2.*nuCold*tr(EPS(uxCold, uyCold, uzCold))*EPS(uxCold,uyCold,uzCold) 
                        + 2.*nuCold*tr(EPS(uxCold, uyCold, uzCold))*EPS(vxDCold,vyDCold,vzDCold)  )  
                        *(-thetaxp*nx - thetayp*ny - thetazp*nz)
                 ) 
     +on(1,2,3,4,5,6,thetax=0,thetay=0,thetaz=0);

varf shapeDerivativeDHot([thetax,thetay,thetaz],[thetaxp,thetayp,thetazp])
     =int2d(ThBox,10,optimize=0)( 
                        (-2.*nuHot*tr(EPS(uxHot, uyHot, uzHot))*EPS(uxHot,uyHot,uzHot) 
                         +2.*nuHot*tr(EPS(uxHot, uyHot, uzHot))*EPS(vxDHot,vyDHot,vzDHot)  )  
                        *(thetaxp*nx + thetayp*ny + thetazp*nz)
                 ) 
     +on(1,2,3,4,5,6,thetax=0,thetay=0,thetaz=0);
     
real[int] diffJ=shapeDerivativeJ(0,Vh,tgv=-1);
real[int] diffVHot=shapeDerivativeVHot(0,Vh,tgv=-1);
real[int] diffDCold=shapeDerivativeDCold(0,Vh,tgv=-1);
real[int] diffDHot=shapeDerivativeDHot(0,Vh,tgv=-1);

real gammaR = 10. * $hmin;

real r = 0.5; real zmin=-5.5; real zmax=5.5; 
func cut=min(max((x+1.)^2 + (y+2.25)^2 - r^2, z-zmin), max((x+1.)^2 + (y-2.25)^2 - r^2,zmax-z) );
Ph cutp1=cut<0;

varf riesz([thetax,thetay,thetaz],[thetaxp,thetayp,thetazp])=int3d(ThBox)(gammaR^2*GRAD(thetax,thetay,thetaz)'*GRAD(thetaxp,thetayp,thetazp)+[thetax,thetay,thetaz]'*[thetaxp,thetayp,thetazp])                                                            
+int3d(ThBox)(tgv*cutp1*[thetax,thetay,thetaz]'*[thetaxp,thetayp,thetazp])
+on(1,2,3,4,5,6,thetax=0,thetay=0,thetaz=0);

Mat A; //Global distributed matrix
{
macro def(i)[i, i#B, i#C] //EOM
macro init(i)[i, i, i] // EOM
createMat(ThBox, A, [P1,P1,P1])
}

A = riesz(Vh, Vh, tgv=-1);
set(A, sparams = "-pc_type hypre"); //"-pc_type lu"

thetaxJ[] = A^-1 * diffJ;
thetaxVHot[] = A^-1 * diffVHot;
thetaxDCold[] = A^-1 * diffDCold;
thetaxDHot[] = A^-1 * diffDHot;

fespace VhGlobal(Th, [P1,P1,P1]);
VhGlobal [thetaxJGlob,thetayJGlob,thetazJGlob], [thetaxDColdGlob,thetayDColdGlob,thetazDColdGlob], [thetaxDHotGlob,thetayDHotGlob,thetazDHotGlob], [thetaxVHotGlob,thetayVHotGlob,thetazVHotGlob];
VhGlobal [thetaxJGlobSum,thetayJGlobSum,thetazJGlobSum], [thetaxDColdGlobSum,thetayDColdGlobSum,thetazDColdGlobSum], [thetaxDHotGlobSum,thetayDHotGlobSum,thetazDHotGlobSum], [thetaxVHotGlobSum,thetayVHotGlobSum,thetazVHotGlobSum];
int[int] subIdx = restrict(Vh, VhGlobal, n2o);

//I want to convert gradJ from local to global!

thetaxJ[] .*= A.D;
thetaxJGlob[](subIdx) = thetaxJ[];
mpiReduce(thetaxJGlob[], thetaxJGlobSum[], processor(0, mpiCommWorld), mpiSUM);

thetaxDCold[] .*= A.D;
thetaxDColdGlob[](subIdx) = thetaxDCold[];
mpiReduce(thetaxDColdGlob[], thetaxDColdGlobSum[], processor(0,mpiCommWorld), mpiSUM);

thetaxDHot[] .*= A.D;
thetaxDHotGlob[](subIdx) = thetaxDHot[];
mpiReduce(thetaxDHotGlob[], thetaxDHotGlobSum[], processor(0, mpiCommWorld), mpiSUM);

thetaxVHot[] .*= A.D;
thetaxVHotGlob[](subIdx) = thetaxVHot[];
mpiReduce(thetaxVHotGlob[], thetaxVHotGlobSum[], processor(0, mpiCommWorld), mpiSUM);

//1st option
real[int] diffJGlob(VhGlobal.ndof), diffDColdGlob(VhGlobal.ndof), diffDHotGlob(VhGlobal.ndof), diffVHotGlob(VhGlobal.ndof);
real[int] diffJGlobSum(VhGlobal.ndof), diffDColdGlobSum(VhGlobal.ndof), diffDHotGlobSum(VhGlobal.ndof), diffVHotGlobSum(VhGlobal.ndof);

diffJ .*= A.D;
diffJGlob(subIdx) = diffJ;
mpiReduce(diffJGlob, diffJGlobSum, processor(0, mpiCommWorld), mpiSUM); 

diffDCold .*= A.D;
diffDColdGlob(subIdx) = diffDCold;
mpiReduce(diffDColdGlob, diffDColdGlobSum, processor(0, mpiCommWorld), mpiSUM);

diffDHot .*= A.D;
diffDHotGlob(subIdx) = diffDHot;
mpiReduce(diffDHotGlob, diffDHotGlobSum, processor(0, mpiCommWorld), mpiSUM);

diffVHot .*= A.D;
diffVHotGlob(subIdx) = diffVHot;
mpiReduce(diffVHotGlob, diffVHotGlobSum, processor(0, mpiCommWorld), mpiSUM);


if (mpirank==0){
    saveArray("$OUTPUT/diffJ.gp",diffJGlobSum);
    saveArray("$OUTPUT/diffVHot.gp",diffVHotGlobSum);
    saveArray("$OUTPUT/diffDCold.gp",diffDColdGlobSum);
    saveArray("$OUTPUT/diffDHot.gp",diffDHotGlobSum);

    saveArray("$OUTPUT/gradJ.gp",thetaxJGlobSum[]);
    saveArray("$OUTPUT/gradVHot.gp",thetaxVHotGlobSum[]);
    saveArray("$OUTPUT/gradDCold.gp",thetaxDColdGlobSum[]);
    saveArray("$OUTPUT/gradDHot.gp",thetaxDHotGlobSum[]);

    //savevtk("$OUTPUT/gradJ.vtu",Th, [-thetaxJGlobSum,-thetayJGlobSum,-thetazJGlobSum]);
    //savevtk("$OUTPUT/gradDCold.vtu",Th, [-thetaxDColdGlobSum,-thetayDColdGlobSum,-thetazDColdGlobSum]);
    //savevtk("$OUTPUT/gradDHot.vtu",Th, [-thetaxDHotGlobSum,-thetayDHotGlobSum,-thetazDHotGlobSum]);
    //savevtk("$OUTPUT/gradVHot.vtu",Th, [-thetaxVHotGlobSum,-thetayVHotGlobSum,-thetazVHotGlobSum]);
}
'''

def sensitivity(mesh):
    FreeFemRunner([preamble, code_sensitivity],config,run_dir=output,run_file='sensitivities.edp',debug=1).execute({'MESH':mesh, 'hmin':hmin},ncpu=N)

    dJ=readFFArray(output+'/diffJ.gp')
    dDc=readFFArray(output+'/diffDCold.gp')
    dDh=readFFArray(output+'/diffDHot.gp')
    dVh=readFFArray(output+'/diffVHot.gp')

    gradJ=readFFArray(output+'/gradJ.gp')
    gradDc=readFFArray(output+'/gradDCold.gp')
    gradDh=readFFArray(output+'/gradDHot.gp')
    gradVh=readFFArray(output+'/gradVHot.gp')

    return (dJ,dDc,dDh,dVh,gradJ,gradDc,gradDh,gradVh)


def shapeUpdate(mesh, xi):
    M = Mesh3D(mesh)
    phi = mshdist(M,ncpu=N) # Initial level set

    # Advection of the level set
    vel = P1Vector3D(M,[xi[::3],xi[1::3],xi[2::3]])
    #phiNew = advect(M,phi,vel,T=1.) #T=1.
    M.save(output+"/M.meshb")
    phi.save(output+"/phi.sol")
    vel.save(output+"/vel.sol")
    os.system(f"advect -dt 1. -nocfl {output}/M.meshb -s {output}/vel.sol -c {output}/phi.sol -o {output}/phiNew.solb -v")
    phiNew = P1Function3D(M,output+"/phiNew.solb")

    ##CHANGE THIS PART, PUT the following barrier 
    # Enforce non optimizable regions
    cutP1 = P1Function3D(M,cut)
    phiNew.sol = np.maximum(phiNew.sol,cutP1.sol)
    barrierP1 = P1Function3D(M,barrier)
    phiNew.sol = np.minimum(phiNew.sol, barrierP1.sol)

    newM=mmg3d(M,hmin,hmax,hgrad,hausd,ls=True,sol=phiNew,params=refermmg,extra_args="-rmc 1e-4",debug=1)
    newM.save(output+'/Th.o.mesh')
    
    d1 = mshdist(newM,ncpu=N)
    d1.save(output+"/d1.sol")
    return output+'/Th.o.mesh'


class HeatExchanger(Optimizable):
    def __init__(self):
        super().__init__()
        self.ucomputed = False
        self.sensitivity_computed = False
        self.obj = None
        self.nconstraints = 1
        self.nineqconstraints = 2

    # Initialization
    def x0(self):
        return output+"/Th0.mesh"

    # Objective function and constraints
    def evalObjective(self, x):
        if not self.ucomputed:
            (J,Dc,Dh,Vh)=computeU(x)
            self.ucomputed = True
            self.obj = (J,Dc,Dh,Vh)
        return self.obj

    def J(self,x):
        (J,Dc,Dh,Vh)=self.evalObjective(x)
        return J

    def H(self,x):
        (J,Dc,Dh,Vh)=self.evalObjective(x)
        return [Dc,Dh]

    def G(self,x):
        (J,Dc,Dh,Vh)=self.evalObjective(x)
        return [Vh]

    # Shape derivatives, sensitivity of objective and constraint
    def evalSensitivities(self,x):
        if not self.sensitivity_computed:
            (dJ,dDc,dDh,dVh,gradJ,gradDc,gradDh,gradVh)=sensitivity(x)
            self.sensitivity_computed = True
            self.sensitivities = (dJ,dDc,dDh,dVh,gradJ,gradDc,gradDh,gradVh)
        return self.sensitivities

    def dJ(self,x):
        (dJ,dDc,dDh,dVh,gradJ,gradDc,gradDh,gradVh) = self.evalSensitivities(x)
        return dJ

    def dH(self,x):
        (dJ,dDc,dDh,dVh,gradJ,gradDc,gradDh,gradVh) = self.evalSensitivities(x)
        return [dDc,dDh]

    def dG(self,x):
        (dJ,dDc,dDh,dVh,gradJ,gradDc,gradDh,gradVh) = self.evalSensitivities(x)
        return [dVh]

    #Gradient and transpose
    def dJT(self,x):
        (dJ,dDc,dDh,dVh,gradJ,gradDc,gradDh,gradVh) = self.evalSensitivities(x)
        return gradJ

    def dHT(self,x):
        (dJ,dDc,dDh,dVh,gradJ,gradDc,gradDh,gradVh) = self.evalSensitivities(x)
        return np.asarray([gradDc,gradDh]).T

    def dGT(self,x):
        (dJ,dDc,dDh,dVh,gradJ,gradDc,gradDh,gradVh) = self.evalSensitivities(x)
        return np.asarray([gradVh]).T

    #Retraction : shape update
    def retract(self, x, dx):
        # Assume that the next computations will be performed on a new mesh
        self.sensitivity_computed = False
        self.ucomputed = False
        return shapeUpdate(x,dx)

    #Accept step : save the current result
    def accept(self,results):
        mesh = results['x'][-1]
        niter = len(results['J'])-1
        itf = format(niter,'04d')
        shutil.copyfile(mesh,output+'/Th_'+itf+'.mesh')
        shutil.copyfile(output+'/d1.sol',output+'/d1.o.sol')
        shutil.copyfile(output+'/THot_64_0000_00.vtu',output+'/THot_'+itf+'.vtu')
        shutil.copyfile(output+'/TCold_64_0000_00.vtu',output+'/TCold_'+itf+'.vtu')
        results['x'][-1] = output+'/Th_'+itf+'.mesh'

        open(output+'/J.txt', 'w').close()
        open(output+'/H.txt', 'w').close()
        open(output+'/G.txt', 'w').close()

        df = pd.DataFrame(data=np.array(results['J']),columns=['J'])
        df.to_csv(output+'/J.txt', header=None, index=None, sep=' ', mode='a')
        df = pd.DataFrame(data=np.array(results['H']),columns=['Dc','Dh'])
        df.to_csv(output+'/H.txt', header=None, index=None, sep=' ', mode='a')
        df = pd.DataFrame(data=np.array(results['G']),columns=['Vh'])
        df.to_csv(output+'/G.txt', header=None, index=None, sep=' ', mode='a')

optSettings = {'dt':hmin,'alphaJ':2.,'alphaC':1., 'maxit':500, 'provide_gradient': True, 'maxtrials':1,
              'itnormalisation': 3}
results=nlspace_solve(HeatExchanger(), optSettings)
