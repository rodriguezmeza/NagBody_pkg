	DOUBLE PRECISION t0
	write(*,*) 'testing time given by the cpu'
	CALL cputime(t0)
	write(*,*) 't0=',t0
	stop
	end
