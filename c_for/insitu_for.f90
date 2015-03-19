program main
use mpi
use iso_c_binding, only: C_NULL_CHAR, C_CHAR, C_INT, C_DOUBLE, C_PTR, C_LOC

implicit none
interface
     !subroutine copyfld (pid, timestep, fldname, startdim1, enddim1, startdim2, enddim2, field) bind(C)
     subroutine cpdata (pid, timestep, fldname, numdims, dimextents, field, com_group) bind(C)
       use, intrinsic :: iso_c_binding
       integer(C_INT), intent(in), value :: pid
       integer(C_INT), intent(in), value :: timestep
       character(len=1), intent(in)      :: fldname
       integer(C_INT), intent(in), value :: numdims
       integer(C_INT), intent(in)        :: dimextents(3)
     !  integer(C_INT), intent(in), value :: startdim1
     !  integer(C_INT), intent(in), value :: enddim1
     !  integer(C_INT), intent(in), value :: startdim2
     !  integer(C_INT), intent(in), value :: enddim2
       real(C_DOUBLE), intent(in)        :: field(4,4)
       integer(C_INT), intent(in) :: com_group
      !type(C_PTR), intent(in), value 		:: fld
     end subroutine cpdata
     !end subroutine copyfld
end interface

integer ierr, rank, numrank
character(8) :: name1 = "Pressure"
character(4) :: extension = ".txt"
character(10) :: format_string
character(20) :: filename
character(10) :: arg
integer counter, time_s
integer extent, i, j
real*8, dimension(:,:), allocatable :: field 
integer, dimension(3) :: dimextents
counter=0
format_string = "(A8,I1,A4)"
name1 = "Pressure"
extension = ".txt"
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numrank, ierr)

if( rank .eq. 0) then
  call getarg(1,arg)
  read(arg, "(i1)") extent
  dimextents(1:2) = extent
  dimextents(3) = 1
end if

call MPI_BCAST(dimextents, 3, MPI_INT, 0, MPI_COMM_WORLD, ierr)

write (filename,format_string) name1, rank, extension
allocate(field(dimextents(1), dimextents(2)))

open(9, file=trim(filename))
do j=1, dimextents(1)
  do i=1, dimextents(2)
    read(9,*) field(i,j)
    counter=counter+1
  end do
end do
!write(*,*) counter
close(9)
!write(*,*) rank
time_s=0;
do time_s=0, 10
  call cpdata (rank, time_s, "T", 2, dimextents(:), field(1:4, 1:4), MPI_COMM_WORLD)
end do

call MPI_FINALIZE(ierr)
end program main
