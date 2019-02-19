program test
   use types_mod
   implicit none
   integer :: ngpt,natoms,i
   real*8 :: dens(3),vnuc(3),q(4),pos(3,4),cell(3,3),dens2(2)
   type(T_Molecule) :: mol1,mol2

   ngpt=3
   dens=(/1,2,3/)
   vnuc=(/4,5,6/)
   natoms=4
   q=(/1,2,3,4/)
   do i=1,4
      pos(1,i)=i*0.1
      pos(2,i)=i*0.2
      pos(3,i)=i*0.3
   enddo
   do i=1,3
      cell(i,i)=5
   enddo

   mol1=T_Molecule(ngpt=ngpt,density=dens)
   mol2=mol1

   call Mol_set_vnuc(mol1,vnuc)
   call Mol_init_ions(mol1,natoms,q,pos)
   call Mol_init_grid(mol1,cell=cell)

   write(*,*) Mol_has_vnuc(mol1),Mol_has_ions(mol1),Mol_has_grid(mol1)
   write(*,*) mol1%vnuc
   write(*,*) Ions_has_charges(mol1%ions)
   write(*,*) mol1%grid%cell
   write(*,*) Mol_has_vnuc(mol2),Mol_has_ions(mol2),Mol_has_grid(mol2)

   call Mol_set_ions(mol2,mol1%ions)
   call Mol_set_grid(mol2,mol1%grid)

   write(*,*) Mol_has_vnuc(mol2),Mol_has_ions(mol2),Mol_has_grid(mol2)

   dens2=(/1,2/)
   call Mol_set_density(mol2,dens2)

end program