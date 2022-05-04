c**********************************************************************
c     matmul.f - matrix - vector multiply, simple self-scheduling version
c************************************************************************
     
      Program Matmult
c#######################################################################
c#
c# This is an MPI example of multiplying a vector times a matrix
c# It demonstrates the use of :
c#
c# * MPI_Init
c# * MPI_Comm_rank
c# * MPI_Comm_size
c# * MPI_Bcast
c# * MPI_Recv
c# * MPI_Send
c# * MPI_Finalize
c# * MPI_Abort
c#
c#######################################################################
    

      include 'mpif.h'

      integer MAX_ROWS, MAX_COLS, rows, cols
      parameter (MAX_ROWS = 1000, MAX_COLS = 1000, MAX_PROCS =32)
      double precision a(MAX_ROWS,MAX_COLS), b(MAX_COLS), c(MAX_COLS)
      double precision buffer(MAX_COLS), ans
      integer procs(MAX_COLS), proc_totals(MAX_PROCS)

      integer myid, master, numprocs, ierr, status(MPI_STATUS_SIZE)
      integer i, j, numsent, numrcvd, sender, job(MAX_ROWS)
      integer rowtype, anstype, donetype

      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      if (numprocs .lt. 2) then
         print *, "Must have at least 2 processes!"
         call MPI_ABORT( MPI_COMM_WORLD, 1 )
         stop
      else if (numprocs .gt. MAX_PROCS) then
         print *, "Must have 32 processes or less."
         call MPI_ABORT( MPI_COMM_WORLD, 1 )
         stop        
      endif
      print *, "Process ", myid, " of ", numprocs, " is alive"

      rowtype  = 1
      anstype  = 2
      donetype = 3

      master   = 0
      rows     = 100
      cols     = 100

      if ( myid .eq. master ) then
c        master initializes and then dispatches
c        initialize a and b
         do 20 i = 1,cols
            b(i) = 1
            do 10 j = 1,rows
               a(i,j) = i
 10         continue
 20      continue

         numsent = 0
         numrcvd = 0
         
c        send b to each other process
         call MPI_BCAST(b, cols, MPI_DOUBLE_PRECISION, master,
     $        MPI_COMM_WORLD, ierr)

c        send a row to each other process
         do 40 i = 1,numprocs-1
            do 30 j = 1,cols
               buffer(j) = a(i,j)
 30         continue
            call MPI_SEND(buffer, cols, MPI_DOUBLE_PRECISION, i,
     $           rowtype, MPI_COMM_WORLD, ierr)
            job(i)  = i
            numsent = numsent+1
 40      continue
         
         do 70 i = 1,rows
            call MPI_RECV(ans, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE,
     $           anstype, MPI_COMM_WORLD, status, ierr)
            sender = status(MPI_SOURCE)
            c(job(sender)) = ans
            procs(job(sender))= sender
            proc_totals(sender+1) =  proc_totals(sender+1) +1

            if (numsent .lt. rows) then
               do 50 j = 1,cols
                  buffer(j) = a(numsent+1,j)
 50            continue
               call MPI_SEND(buffer, cols, MPI_DOUBLE_PRECISION, sender,
     $              rowtype, MPI_COMM_WORLD, ierr)
               job(sender) = numsent+1
               numsent     = numsent+1
            else
            call MPI_SEND(1, 1, MPI_INTEGER, sender, donetype,
     $           MPI_COMM_WORLD, ierr)
            endif
 70      continue
         
c        print out the answer
        
         do 80 i = 1,cols
	    write(6,809) i,c(i),procs(i)
 809	    format('c(',i3,') =',f8.2,' computed by proc #',i3)
 80      continue
         do 81  i=1,numprocs
            write(6,810) i-1,proc_totals(i)
 810	    format('Total answers computed by processor #',i2,' were ',i3)
 81        continue

      else
c        compute nodes receive b, then compute dot products until done message
         call MPI_BCAST(b, cols, MPI_DOUBLE_PRECISION, master,
     $        MPI_COMM_WORLD, ierr)
 90      call MPI_RECV(buffer, cols, MPI_DOUBLE_PRECISION, master,
     $        MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
         if (status(MPI_TAG) .eq. donetype) then
            go to 200
         else
            ans = 0.0
            do 100 i = 1,cols
               ans = ans+buffer(i)*b(i)
 100        continue
            call MPI_SEND(ans, 1, MPI_DOUBLE_PRECISION, master, anstype,
     $           MPI_COMM_WORLD, ierr)
            go to 90
         endif
      endif

 200  call MPI_FINALIZE(ierr)
      stop
      end




