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

   call mol1%set_vnuc(vnuc)
   call mol1%init_ions(natoms,q,pos)
   call mol1%init_grid(cell=cell)

   write(*,*) mol1%has_vnuc(),mol1%has_ions(),mol1%has_grid()
   write(*,*) mol1%vnuc
   write(*,*) mol1%ions%has_charges()
   write(*,*) mol1%grid%cell
   write(*,*) mol2%has_vnuc(),mol2%has_ions(),mol2%has_grid()

   call mol2%set_ions(mol1%ions)
   call mol2%set_grid(mol1%grid)

   write(*,*) mol2%has_vnuc(),mol2%has_ions(),mol2%has_grid()

   dens2=(/1,2/)
   call mol2%set_density(dens2)

end program
