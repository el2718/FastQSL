!contact: Chen, Jun; el2718@mail.ustc.edu.cn
!gfortran -O3 sudoku.f90 -o sudoku.x

module share
implicit none
integer(8)::n_solved, solved_max
integer::eliminated_max, m_shift(9,9)
character(len=127)::filename
logical::outfile, verbose, brute_force
end module share


program main
use share
implicit none
integer::sudoku(9,9), sudoku_orig(9,9), i, j, narg
real::x(9,9)
character(len=127):: puzzle, line_str, arg
logical:: exist_flag
!--------------------------------------------
! default setting
brute_force=.false.
outfile=.false.
solved_max=2
eliminated_max=0

narg=iargc()

do i=1, narg
	call get_command_argument(i, arg)
	if (arg .eq. '-h') then
		print*, 'usage: ./sudoku.x [puzzle] [-h] [-b] [-o] [-s solved_max] [-e eliminated_max]'
		print*, '' 
		print*, 'puzzle: a text file with 81 integer elements of a puzzle, a 0 or a . represents an empty cell'
		print*, '' 
		print*, '-h, find help' 
		print*, '' 
		print*, '-b, use brute force to solve the puzzle; just try and check the consistency with back tracking.  Otherwise logical strategies are used by default, see https://www.sudokuwiki.org/Strategy_Families'
		print*, '' 
		print*, '-o, export results to text files. A puzzle with 81 non-empty elements will not be exported'
		print*, '' 
		print*, '-s solved_max; the maximum number of solutions that the program could give, in case a sudoku may have multiple solutions.'
		print*, '   --If set to be larger than 1, non-repeative solutions will be presented'
		print*, '   --If set to a number more than the count of all solutions, the count will be reported'
		print*, '   --The default value of solved_max is 2, to check the uniqueness of the solution'
		print*, '   --If set to 1, will not check the uniqueness, but will be faster'
		print*, '' 
		print*, '-e eliminated_max; If a puzzle has a unique solution, eliminating some non-empty elements from the puzzle could still give the same solution. This value is the maximum number of non-empty elements can be eliminated for the same unique solution. If then no any more non-empty element can be eliminated, the actual eliminated number could be smaller than this value'
		print*, '   --The default value of eliminated_max is 0. If eliminated_max is set to > 0, the eliminated puzzle has the same unique solution will be present. The sequence of elimination is random. This can create a puzzle from a full filled sudodu, e.g. from a solution of another puzzle'
		return
	endif
	if (trim(arg) .eq. '-b') brute_force=.true.
	if (trim(arg) .eq. '-o') outfile=.true.

	if (trim(arg) .eq. '-s' .and. i .ne. narg) then
		call get_command_argument(i+1, arg)
		read(arg,'(i19)') solved_max
	endif
	if (trim(arg) .eq. '-e' .and. i .ne. narg) then
		call get_command_argument(i+1, arg)
		read(arg, '(i9)') eliminated_max
	endif
enddo

call get_command_argument(1, puzzle)

if (puzzle .eq. "" .or.         &
	trim(puzzle) .eq. "-h" .or. &
    trim(puzzle) .eq. "-b" .or. &
	trim(puzzle) .eq. "-o" .or. &
	trim(puzzle) .eq. "-s" .or. &
	trim(puzzle) .eq. "-e") then
	sudoku=0
	if (outfile) filename='All_zero'
else
	inquire(file = trim(puzzle), exist=exist_flag)

	if (exist_flag) then
		if (outfile) then
			do i=len(puzzle),1,-1
				if (puzzle(i:i) .eq. ".") then
					filename=puzzle(1:i-1)
					exit
				endif
			enddo
		endif

		open(unit=8, file=trim(puzzle), status='old')
		do j=1, 9
			read(8, "(A)") line_str
			do i=1,len(line_str)
				if (line_str(i:i) .eq. '.') line_str(i:i)='0'
			enddo
			read(line_str,*) sudoku(:,j)
		enddo
		close(8)
	else
		print*, trim(puzzle)//' not exist'
		return
	endif
endif

call random_seed()
call random_number(x)
m_shift=floor(x*9)
!--------------------------------------------
write(*,"('  == input ==')")
call print_sudoku(sudoku)
sudoku_orig=sudoku
verbose=.true.
call resolve(sudoku)	

if (n_solved .eq. 1 .and. solved_max .ge. 2 .and. eliminated_max .gt. 0) call eliminate(sudoku_orig)

end program main


subroutine resolve(sudoku)
use share
implicit none
integer::sudoku(9,9)
logical::candidate(9,9,9), bug_flag
!--------------------------------------------
call check_sudoku(sudoku, bug_flag)
if (bug_flag) then
	write(*,"('problematic input')")
	return
endif
!--------------------------------------------
n_solved=0
if (any(sudoku .eq. 0)) then
	if (brute_force) then
		!brute force
		call try_sudoku(1, sudoku)
	else
		!logical strategies
		call initialize_candidate(sudoku, candidate)
		call process(candidate, bug_flag)
		if (.not. bug_flag)	call try_candidate(candidate, bug_flag)
	endif
else
	n_solved=1
endif
!-------------------------------------------
if (n_solved .lt. solved_max .and. verbose) then
	select case (n_solved)
	case(0)
		write(*,"('  == it do not have any solution ==')")
	case(1)	
		write(*,"('  == it has a unique solution ==')")
	case default
		write(*,"('  == it has ', i0,' solutions ==')") n_solved
	end select
endif
!--------------------------------------------
end subroutine resolve


subroutine process(candidate, bug_flag)
implicit none
integer::no_update_times, n_candidate, n_candidate0
logical::candidate(9,9,9), bug_flag
!--------------------------------------------
no_update_times=0
n_candidate0=count(candidate)
do while(count(candidate)>81)
	call strategy(no_update_times+1, candidate, bug_flag)	
	if (bug_flag) return
	
	n_candidate=count(candidate)	
	if(n_candidate0==n_candidate) then
		no_update_times=no_update_times+1
		if (no_update_times==3) return
	else
		n_candidate0=n_candidate
		no_update_times=0
	endif
enddo
end subroutine process


subroutine strategy(n, candidate, bug_flag)
!https://www.sudokuwiki.org/strategy_families
!n=1, basic
!n=2, naked pairs, hidden pairs
!n=3, naked triples, hidden triples
implicit none
integer::k, n, m, way, group(9), i_group, c9n, i_combination
integer::i_series(9), j_series(9), indgen(9)=(/1,2,3,4,5,6,7,8,9/)
logical::candidate(9,9,9), candidates(9), positions(9), bug_flag
!--------------------------------------------
call combination_number(9, n, c9n)
do i_combination=1,c9n	
	call combination_group(n, group, i_combination)
do way=1,3
do k=1,9
!--------------------------------------------
!if n positions have only n candidates
	call ij_series(k, way, group, i_series, j_series)
	candidates=candidate(:, i_series(1), j_series(1))	
	do i_group=2,n		
		candidates=candidates .or. candidate(:,i_series(i_group),j_series(i_group))	
	enddo
	
	bug_flag=count(candidates) .lt. n
	if(bug_flag) return	
	
	if (count(candidates) .eq. n) then
		do i_group=n+1,9
			where(candidates) candidate(:,i_series(i_group),j_series(i_group))=.false.
		enddo
	endif	
!--------------------------------------------
!if n candidates found in only n positions
	call ij_series(k, way, indgen, i_series, j_series)
	forall(m=1:9) positions(m)=any(candidate(group(1:n), i_series(m), j_series(m)))
	
	bug_flag=count(positions) .lt. n
	if(bug_flag) return
	
	if (count(positions) .eq. n) &
	forall(m=1:9, positions(m)) candidate(group(n+1:9), i_series(m), j_series(m))=.false.
!--------------------------------------------
enddo
enddo
enddo
call check_candidate(candidate, bug_flag)
end subroutine strategy


recursive subroutine try_candidate(candidate, bug_flag)
use share
implicit none
integer::sudoku(9,9), candidate_first, m, k, i, j, n_solved0, n_guess
logical::candidate(9,9,9), candidate_try(9,9,9), bug_flag
!--------------------------------------------
100 continue

if (count(candidate) .le. 81) then
	call check_candidate(candidate, bug_flag)
	if (bug_flag) return
	forall(i=1:9,j=1:9)	sudoku(i,j)=findloc(candidate(:,i,j),.true.,1)
	n_solved=n_solved+1
	if (verbose) call print_sudoku(sudoku)
	return
endif
!--------------------------------------------
do n_guess=2,9
do j=1,9
do i=1,9
if (count(candidate(:,i,j)) .eq. n_guess) then	
	do m=0,8
		if (candidate(mod(m+m_shift(i,j),9)+1,i,j)) then
			candidate_first=mod(m+m_shift(i,j),9)+1
			exit
		endif	
	enddo
	candidate_try=candidate
	candidate_try(:,i,j)=.false.
	candidate_try(candidate_first,i,j)=.true.
	call process(candidate_try, bug_flag)
	
	n_solved0=n_solved
	if (.not. bug_flag)	call try_candidate(candidate_try, bug_flag)
	
	if (n_solved .eq. solved_max) then
		candidate=candidate_try
		return
	endif
		
	if (bug_flag .or. n_solved>n_solved0) then
		candidate(candidate_first,i,j)=.false.
		call process(candidate, bug_flag)
		
		if (bug_flag) then
			return
		else
			goto 100
		endif
	endif
endif 
enddo
enddo
enddo
end subroutine try_candidate


recursive subroutine try_sudoku(ij, sudoku)
use share
implicit none
integer::sudoku(9,9), sudoku_try(9,9), ij, i, j, m, m0
logical::ij_mark(0:9)
integer::neighbor_i(2), neighbor_j(2)
!--------------------------------------------
if (ij .eq. 82) then
	n_solved=n_solved+1
	if (verbose) call print_sudoku(sudoku)
	return
endif
!--------------------------------------------
j=(ij-1)/9+1
i=mod(ij-1,9)+1
if(sudoku(i,j)==0)then
	sudoku_try=sudoku

	ij_mark=.false.
	ij_mark(sudoku(:,j))=.true.
	ij_mark(sudoku(i,:))=.true.
	call i_neighbor(i,neighbor_i)
	call i_neighbor(j,neighbor_j)
	ij_mark(sudoku(neighbor_i,neighbor_j(1)))=.true.
	ij_mark(sudoku(neighbor_i,neighbor_j(2)))=.true.

	do m0=0,8
	
		m= mod(m0+m_shift(i,j),9)+1
		if (ij_mark(m)) cycle

		sudoku_try(i,j)=m
		call try_sudoku(ij+1, sudoku_try)
		
		if (n_solved .eq. solved_max) then
			sudoku=sudoku_try
			return
		endif
	enddo
else
	call try_sudoku(ij+1, sudoku)
endif
end subroutine try_sudoku


subroutine i_neighbor(i, neighbor)
implicit none
integer::i, neighbor(2)
select case(i)
case(1)
	neighbor=[2,3]
case(2)
	neighbor=[1,3]
case(3)
	neighbor=[1,2]
case(4)
	neighbor=[5,6]
case(5)
	neighbor=[4,6]
case(6)
	neighbor=[4,5]
case(7)
	neighbor=[8,9]
case(8)
	neighbor=[7,9]
case(9)
	neighbor=[7,8]
end select
end subroutine i_neighbor


subroutine eliminate(sudoku)
use share
implicit none
integer::sudoku(9,9), sudoku_try(9,9), &
ij0, ij, ij_shift, i, j, n_eliminated
logical::tested_mark(9,9), bug_flag 
real::x
character(len=2)::n_eliminated_str
character(len=1)::plural_suffix
character(len=127)::file_saved
!--------------------------------------------
n_eliminated=0
tested_mark=.false.
solved_max=2
verbose=.false.
!--------------------------------------------
200 continue

call random_seed()
call random_number(x)
ij_shift=floor(x*81)
!--------------------------------------------
do ij0=0,80
	ij=mod(ij0+ij_shift,81)
	j=ij/9+1
	i=mod(ij,9)+1

	sudoku_try=sudoku	
	if (sudoku(i,j) .gt. 0 .and. .not. tested_mark(i,j)) then	
		tested_mark(i,j)=.true.	
		sudoku_try(i,j)=0
		call resolve(sudoku_try)
		if (n_solved .eq. 1) then
			sudoku(i,j)=0
			n_eliminated=n_eliminated+1
			if(n_eliminated .eq. eliminated_max) exit
			goto 200
		endif
	endif
enddo
!--------------------------------------------
if (n_eliminated .eq. 0) then
	print*, ' == no element can be eliminated =='
	return
endif
!--------------------------------------------
write(n_eliminated_str,"(i0)") n_eliminated
if (n_eliminated .eq. 1) then 
	plural_suffix=""
else
	plural_suffix="s"
endif

print*, ' == eliminate '//trim(n_eliminated_str)&
//' element'//trim(plural_suffix)//' =='

n_solved=0
call print_sudoku(sudoku)
!--------------------------------------------
if (outfile) then
	file_saved=trim(filename)//'_eliminate'//trim(n_eliminated_str)//'.txt'
	print*, trim(file_saved)//' is saved'
	open(unit=8, file=trim(file_saved), status='replace')
	write(8, '(9i2)') ((sudoku(i,j),i=1,9),j=1,9)
	close(8)
endif
end subroutine eliminate


subroutine initialize_candidate(sudoku, candidate)
implicit none
integer::sudoku(9,9), i, j
logical::candidate(9,9,9)
!--------------------------------------------
candidate=.true.

do j=1,9
do i=1,9
	if (sudoku(i,j)>0) then
		candidate(:,i,j)=.false.
		candidate(sudoku(i,j),i,j)=.true.
	endif
enddo
enddo
end subroutine initialize_candidate


subroutine check_candidate(candidate, bug_flag)
implicit none
integer::k, m, i0, j0
logical::candidate(9,9,9), bug_flag
!--------------------------------------------
do k=1,9
	j0=(k-1)/3*3+1
	i0=mod(k-1,3)*3+1
do m=1,9
	!check no any candidate in a cell
	bug_flag= .not. (any(candidate(:,m,k)) .and. & 
	!check the candidate m lacked in a unit
	any(candidate(m,k,:)) .and. any(candidate(m,:,k)) .and. any(candidate(m,i0:i0+2,j0:j0+2)))
	if(bug_flag) return
enddo
enddo
end subroutine check_candidate


subroutine check_sudoku(sudoku, bug_flag)
implicit none
integer::sudoku(9,9), k, m, i0, j0
logical::bug_flag
!--------------------------------------------
do k=1,9
	j0=(k-1)/3*3+1
	i0=mod(k-1,3)*3+1
do m=1,9
	bug_flag=count(sudoku(k,:) .eq. m)>1 .or. & 
	         count(sudoku(:,k) .eq. m)>1 .or. & 
	         count(sudoku(i0:i0+2,j0:j0+2) .eq. m)>1
	if(bug_flag) return
enddo
enddo
end subroutine check_sudoku


subroutine print_sudoku(sudoku)
use share
implicit none
integer:: sudoku(9,9), i, j
character(len=30)::row_str, divide_str
character(len=19)::n_solved_str
character(len=127)::file_saved
!--------------------------------------------
if (n_solved .gt. 0) then
	write(n_solved_str,"(i0)") n_solved
	if(outfile) then
		file_saved=trim(filename)//'_solution'//trim(n_solved_str)//'.txt'
		print*, ' == solution '//trim(n_solved_str)//' == '//trim(file_saved)//' is saved'
		open(unit=8, file=trim(file_saved), status='replace')
		write(8, '(9i2)') ((sudoku(i,j),i=1,9),j=1,9)
		close(8)
	else
		print*, ' == solution '//trim(n_solved_str)//' == '
	endif
endif

divide_str='       ------+-------+------'
do j=1,9
	write(row_str,"(6x, 3i2,' |',3i2,' |',3i2)") (sudoku(i,j),i=1,9)
	do i=8,30
		if (row_str(i:i) .eq. '0') row_str(i:i)='.'
	enddo
	write(*,*) row_str	
	if (j==3 .or. j==6) write(*,*) divide_str
enddo
end subroutine print_sudoku


subroutine ij_series(k, way, group, i_series, j_series)
integer::k, way, group(9), i_series(9), j_series(9)
!--------------------------------------------
select case(way)
case(1)
	i_series=k; j_series=group
case(2)	
	i_series=group; j_series=k
case(3)
	j_series=((k-1)/3)*3+1+(group-1)/3
	i_series=(mod(k-1,3))*3+1+mod(group-1,3)
end select
end subroutine ij_series


subroutine combination_number(m, n, cmn)
implicit none
!integer::i
integer::m, n, cmn
!--------------------------------------------
!m=9
if (n .eq. 1 .or. n .eq. 8) cmn=9
if (n .eq. 2 .or. n .eq. 7) cmn=36
if (n .eq. 3 .or. n .eq. 6) cmn=84
if (n .eq. 4 .or. n .eq. 5) cmn=126
! cmn=m-n+1
!do i=m-n+2,m
!	cmn=cmn*i
!enddo
!do i=2,n
!	cmn=cmn/i
!enddo
end subroutine combination_number


subroutine combination_group(n, group, i)
implicit none
integer::i, j, k, n, group(9), step
!--------------------------------------------
do k=1,9
	group(k)=k
enddo

do step=2,i
	call group_plus1(group, n, n)
enddo

do k=n+1,9
do j=1,k-1
	if (all(j .ne. group(1:k-1))) then 
		group(k)=j
		exit
	endif
enddo
enddo
end subroutine combination_group


recursive subroutine group_plus1(group, k, n)
implicit none
integer::group(9), k, n
!--------------------------------------------	
if (group(k) .lt. 9-n+k) then
	group(k)=group(k)+1
else
	if (k .eq. 1) return
	call group_plus1(group, k-1, n)
	group(k)=group(k-1)+1
endif
end subroutine group_plus1
