Requirements:
- petsc version >= 3.11.3 configured with SupoerLU_DIST : https://petsc.org/release/manualpages/Mat/MATSOLVERSUPERLU_DIST/
- FreeFem++ (configured with PETSc): https://doc.freefem.org/introduction/installation.html#compilation-on-ubuntu
- python >= 3.8, with pymedit, pyfreefem, nullspace_optimizer: https://gitlab.com/florian.feppon/pyfreefem, https://gitlab.com/florian.feppon/pymedit, https://gitlab.com/florian.feppon/null-space-optimizer
- ISCDtoolbox (medit, advection, mshdist): https://github.com/iscdtoolbox](https://github.com/ISCDtoolbox/Mshdist, https://github.com/ISCDtoolbox/Advection, https://github.com/ISCDtoolbox/Medit

  Once you have everything installed, you can run the test cases by doing:
  python HETube.py

  Note that the simulations were run on a cluster using 64 processors, you can change the number of processors in the line "N=64"
  and then modify the lines

  shutil.copyfile(output+'/THot_64_0000_00.vtu',output+'/THot_'+itf+'.vtu')
  shutil.copyfile(output+'/TCold_64_0000_00.vtu',output+'/TCold_'+itf+'.vtu')

  replacing 64 by the new number of lines and the 00.vtu by the new digit, for example if you use N=8, then it should read:

  shutil.copyfile(output+'/THot_8_0000_0.vtu',output+'/THot_'+itf+'.vtu')
  shutil.copyfile(output+'/TCold_8_0000_0.vtu',output+'/TCold_'+itf+'.vtu')
