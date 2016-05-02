  PROGRAM vtu2FEHM
!-------------------------------------------------------------------------------
! "vtu2FEHM v.1.0" 
! Copyright 2016 (RSE SpA)
! "vtu2FEHM v.1.0" authors and email contact are provided on 
! the documentation file.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
! This subroutine reads an unstructured grid (.vtu) of a finite element mesh and converts it into a grid input file for FEHM (.grid), 
!    together with the corresponding geological zone list (portion of the file .dat)
! Note: the mesh elements should be all of the same type (element%kind)
! 1) Decoding the unstructured grid (.vtu)
! 2) Writing the grid input file for FEHM (.grid)
! 3) Writing the geological zone list (portion of the file .dat)
! 4) Writing the unstructured .vtu file as a check
! Author: Amicarelli A. (14 Nov 2011)
!
! Declarations 
  type node_type
     real :: coor(3)
     integer :: material
  end type node_type
  type element_type
     integer :: offset,kind,id,material
     integer :: node_id(8)   
  end type element_type
  integer :: n_materials,i,j,k,n_nodes,n_elements,element_nodes,ia,ia_pre
  character(len=80) :: line,char1,char2,vtu_file_name
  integer, dimension(:), allocatable :: material_elem,material_nodes
  character(len=4) :: fmt_char1,fmt_char2
  type(node_type), dimension(:), allocatable :: node
  type(element_type), dimension(:), allocatable :: element
  character(len=80) :: vtu_file_name_test
!
! Initializations
  print *,"!!! Start program Mvtu2FEHM"  
  call getarg (1,vtu_file_name)
  ia_pre=0
  n_materials=0
  char1(1:80)=" "
  char2(1:80)=" "
  vtu_file_name_test = "test.vtu"
!
! 1) Start decoding the unstructured grid (.vtu)
  print *,"!!! 1) Decoding the unstructured grid (.vtu): start"  
  open (1,file=vtu_file_name)
!
! General parameters
  read (1,'(2/)')  
  read (1,'(a)') line
  j=1
  k=0
  do i=1,len(line)
     ia=iachar(line(i:i))
     if ((ia>=48).and.(ia<=57)) then
        if ((ia_pre>=48).and.(ia_pre<=57)) then
           j=j+1
           else
              k=k+1
              j=1
        endif
        if (k==1) then
           char1(j:j)=achar(ia)
           else   
              char2(j:j)=achar(ia)
        endif
     endif
     ia_pre=ia
  end do
  fmt_char1='(i )'
  fmt_char2='(i )'
  write(fmt_char1(3:3),'(i1)') len_trim(char1)
  write(fmt_char2(3:3),'(i1)') len_trim(char2)
  read (char1,fmt_char1) n_nodes 
  read (char2,fmt_char2) n_elements 
  allocate (node(n_nodes))
  allocate (element(n_elements))
  print *,"        The mesh is composed by ", n_nodes, " and ", n_elements, " elements"
!
! Reading the nodes 
  read (1,'(/)')
  do i=1,n_nodes,12
     read(1,'(1x,12(3(f12.4),2x))') (node(j)%coor(:),j=i,i+11)
  end do
  read (1,'(2/)')
! Reading the connectivity structure        
  read (1,'()')         
  do i=1,(n_elements-7),7
     read(1,'(56(1x,i7))') (element(j)%node_id(:),j=i,i+6)
  end do
  do k=i,n_elements
     read(1,'(8(1x,i7))',advance='no') element(k)%node_id(:)
  end do
  read (1,'(/)') 
! Reading the offsets (position of the last element node in the .vtu list)
  read (1,'(a)')   
  do i=1,(n_elements-56),56
     read(1,'(56(1x,i7))') (element(j)%offset,j=i,i+55)
  end do
  do k=i,n_elements
     read(1,'(1x,i7)',advance='no') element(k)%offset
  end do
  read (1,'(/)')
! Reading the element type
  read (1,'()')    
  do i=1,(n_elements-96),96
     read(1,'(96(1x,i2))') (element(j)%kind,j=i,i+95)
  end do
  do k=i,n_elements
     read(1,'(1x,i2)',advance='no') element(k)%kind
  end do
  read (1,'(3/)')
! Reading the mesh element indeces
  read (1,'()')
  do i=1,(n_elements-56),56
     read(1,'(56(1x,i7))') (element(j)%id,j=i,i+55)
  end do
  do k=i,n_elements
     read(1,'(1x,i7)',advance='no') element(k)%id
  end do
  read (1,'(/)')
! Reading the materials
  read (1,'()')
  do i=1,(n_elements-96),96
     read(1,'(96(1x,i2))') (element(j)%material,j=i,i+95)
  end do
  do k=i,n_elements
     read(1,'(1x,i2)',advance='no') element(k)%material
  end do
  do i=1,n_elements
     n_materials=max(n_materials,element(i)%material)
  end do
  allocate (material_elem(n_materials))
  do i=1,n_materials
     material_elem(i)=0
  end do
! Nota: +1 perhè il primo ID dei vtu è zero (convenzione Phyton)
  do i=1,n_elements
     do j=1,8
        element(i)%node_id(j)=element(i)%node_id(j)+1
     end do
  end do
  close(1)
!
  print *,"       Decoding the unstructured grid (.vtu): end" 
! End decoding the unstructured grid (.vtu)
!
! 2) Writing the grid input file for FEHM (.grid)
  print *,"!!! 2) Writing the grid input file for FEHM (.grid): start"  
  open (2,file='FEHM_mesh.grid')
  write (2,'(a)') "coor   n/a     "
  write (2,'(12i)') n_nodes
  j=1
  do i=1,n_nodes
  !AA!!!test
  node(i)%coor(3)=node(i)%coor(3)+2000
  !
  write (2,'(i10,3(f15.5))') i,node(i)%coor(:)
  end do
  write (2,'(a)') "         0        0.00000        0.00000        0.00000" 
  select case (element(1)%kind)
     case (12)
     element_nodes=8
  end select
  write (2,'(a)') "elem trad"
  write (2,'(2(8i))',advance='NO') element_nodes,n_elements
  write (2,'(a)') "   0"
  do i=1,n_elements
  !AA!!!test
!     write (2,'(9(i8))') (i,element(i)%node_id(:))
     write (2,'(9(i8))') i,element(i)%node_id(5),element(i)%node_id(6),element(i)%node_id(7),element(i)%node_id(8), &
                           element(i)%node_id(1),element(i)%node_id(2),element(i)%node_id(3),element(i)%node_id(4)
     material_elem(element(i)%material)=material_elem(element(i)%material)+1
  end do
  write (2,'(a)') "       0       0       0       0       0       0       0       0       0"
  write (2,'(a)') "stop"
  close (2)
  print *,"       Writing the grid input file for FEHM (.grid): end" 
! Writing the grid input file for FEHM (.grid): end
!
! 3) Writing the geological zone list (portion of the file .dat)
  print *,"!!! 3) Writing the geological zone list (portion of the file .dat): start" 
! Each node needs to be given a material
  do i=1,n_nodes
     node(i)%material = 0
     do j=1,n_elements
        do k=1,8 
           if (element(j)%node_id(k) == i) node(i)%material = max(node(i)%material,element(j)%material)   
        end do
     end do
  end do
  allocate (material_nodes(n_materials))
  do i=1,n_materials
     material_nodes(i)=0
  end do
  do i=1,n_nodes
     material_nodes(node(i)%material)=material_nodes(node(i)%material)+1
  end do
  open (3,file='FEHM_zones.dat')
  do i=1,n_materials
     write (3,'(i8)') i
     write (3,'(a)') "nnum"
     write (3,'(i8)') material_nodes(i)
     k=1
     do j=1,n_nodes
        if (node(j)%material==i) then
           k=k+1
           if (MOD(k,56)/=0) then
              write (3,'(i8)',ADVANCE='NO') j
              else
                 write (3,'(i8)') j
           endif
        endif
     end do
     write(3,'()')
  end do
  close (3)
  print *,"       Writing the geological zone list (portion of the file .dat): end" 
! Writing the geological zone list (portion of the file .dat): end
!
! 4) Writing the unstructured .vtu file as a check 
  print *,"!!! 4) Writing the unstructured .vtu file as a check: start  "  
  open (4,file=vtu_file_name_test)
!
! General parameters
  write (4,'(a)') '<?xml version="1.0"?>'
  write (4,'(a)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
  write (4,'(a)') '<UnstructuredGrid>'
  write (4,'(a)') line
  write (4,'(a)') '    <Points>'
  write (4,'(a)') '      <DataArray type="Float32" NumberOfComponents="3" format="ascii" >'
! writing the nodes 
  do i=1,n_nodes,12
     write(4,'(1x,12(3(f12.4),2x))') (node(j)%coor(:),j=i,i+11)
  end do
  write (4,'(a)') '      </DataArray>'
! Writing the connectivity structure        
  write (4,'(a)') '    </Points>'
  write (4,'(a)') '    <Cells>'
  write (4,'(a)') '      <DataArray type="Int32" Name="connectivity" format="ascii" >'
  do i=1,(n_elements-7),7
     write(4,'(56(1x,i7))') (element(j)%node_id(:)-1,j=i,i+6)
  end do
  do k=i,n_elements
     write(4,'(8(1x,i7))',advance='no') element(k)%node_id(:)-1
  end do
  write (4,'()')
  write (4,'(a)') '      </DataArray>'  
! Writing the offsets (position of the last element node in the .vtu list)
  write (4,'(a)') '      <DataArray type="Int32" Name="offsets" format="ascii" >'
  do i=1,(n_elements-56),56
     write(4,'(56(1x,i7))') (element(j)%offset,j=i,i+55)
  end do
  do k=i,n_elements
     write(4,'(1x,i7)',advance='no') element(k)%offset
  end do
  write (4,'()')
  write (4,'(a)') '      </DataArray>'
! Writing the element type
  write (4,'(a)') '      <DataArray type="UInt8" Name="types" format="ascii" >'
  do i=1,(n_elements-96),96
     write(4,'(96(1x,i2))') (element(j)%kind,j=i,i+95)
  end do
  do k=i,n_elements
     write(4,'(1x,i2)',advance='no') element(k)%kind
  end do
  write (4,'()')
  write (4,'(a)') '      </DataArray>'
! Writing the mesh element indeces
  write (4,'(a)') '    </Cells>'
  write (4,'(a)') '    <CellData>'
  write (4,'(a)') '      <DataArray type="Int32" Name="MeshIndex" format="ascii" >'
  do i=1,(n_elements-56),56
     write(4,'(56(1x,i7))') (element(j)%id,j=i,i+55)
  end do
  do k=i,n_elements
     write(4,'(1x,i7)',advance='no') element(k)%id
  end do
  write (4,'()')
  write (4,'(a)') '      </DataArray>'
! Writing the materials
  write (4,'(a)') '      <DataArray type="Int32" Name="Materials" format="ascii" >'
  do i=1,(n_elements-96),96
     write(4,'(96(1x,i2))') (element(j)%material,j=i,i+95)
  end do
  do k=i,n_elements
     write(4,'(1x,i2)',advance='no') element(k)%material
  end do
  write (4,'()')
  write (4,'(a)') '      </DataArray>'
  write (4,'(a)') '    </CellData>'
  write (4,'(a)') '    <PointData>'
  write (4,'(a)') '    </PointData>'
  write (4,'(a)') '  </Piece>'
  write (4,'(a)') '</UnstructuredGrid>'
  write (4,'(a)') '</VTKFile>' 
  close(4)
!
  print *,"       Writing the unstructured .vtu file as a check : end" 
!
  print *,"!!! End program vtu2FEHM" 
  deallocate(node)
  deallocate(material_elem) 
  deallocate(material_nodes) 
  deallocate(element)
  END PROGRAM vtu2FEHM
!





