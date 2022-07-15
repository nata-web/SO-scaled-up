module sofort
  IMPLICIT NONE
CONTAINS
  SUBROUTINE doseed(seed)
    INTEGER*4, INTENT(IN) :: seed
    INTEGER :: size
    INTEGER, allocatable :: seedArray(:)
    call random_seed(size=size)
    allocate(seedArray(size))
    seedArray(:) = seed
    call random_seed(put=seedArray)
	deallocate(seedArray)
  end SUBROUTINE doseed
  
  SUBROUTINE run(w, worig, energies, steps, resets, n, dolearn, eta, speed)
    IMPLICIT NONE

    INTEGER*8, INTENT(IN) :: steps, resets, n, dolearn, eta, speed
    INTEGER*8, DIMENSION(n,n), INTENT(INOUT) :: w
    real*8, DIMENSION(n,n), INTENT(IN) :: worig
    real*8, DIMENSION(steps,resets), INTENT(inOUT) :: energies

    INTEGER*8 :: i
	
    do i=1, resets
      call learn(w, worig, energies(:,i), steps, n, dolearn, eta, speed)
      ! if (mod(i,(resets/10))==0) then
        ! write (*,*) i
      ! end if
    end do
	
	if (speed==0) then
	! if ((speed==0) .and. (dolearn==1)) then
      do i=1, n
        ! Since we only update the upper triangle in learn(), 
        ! we need to copy the other half to make it symmetric again
        w(i+1:n,i) = w(i,i+1:n)
      end do
    end if
	
  end SUBROUTINE run
  
  SUBROUTINE learn(w, worig, energies, steps, n, dolearn, eta, speed)
    IMPLICIT NONE

    INTEGER*8, INTENT(IN) :: steps, n, dolearn, eta, speed
    INTEGER*8, DIMENSION(n,n), INTENT(INOUT) :: w
    real*8, DIMENSION(n,n), INTENT(IN) :: worig
    real*8, DIMENSION(steps), INTENT(inOUT) :: energies

    INTEGER*8 :: t, i, j, idx, delta, activation
    real*8 :: r
    INTEGER*1 :: oldState, newState
    INTEGER*1, DIMENSION(n) :: state
    INTEGER*4, DIMENSION(n) :: idx2t		! Map between idx and the time state(i) was last changed
    INTEGER*4, DIMENSION(steps) :: t2Idx 	! The indices for which at time t the state was changed
    INTEGER*1, DIMENSION(steps) :: t2State	! Time history of all states after they were changed
    INTEGER*8, DIMENSION(n,n) :: dw
	
	dw(:,:) = 0
	
    idx2t(:) = 1 ! ('zeros' in Python)
	
	! Generate a random state
    do i=1,n
      call random_number(r) 	! generate random number from the uniform distribution over the range 0<=x<1
      state(i) = floor(2*r)*2-1 ! convert to -1 or 1
    end do
	
    do t=1,steps
      call random_number(r)
      idx = floor(n*r)+1 	! choose an idx randomly from 1 to n
      oldState = state(idx) ! save the state of that idx
      if (speed==1) then
        if (dolearn==1) then
          w(:,idx) = w(:,idx) + dw(:,idx)*(t-idx2t(idx))
          do i=idx2t(idx)+1, t-1
            newState = state(idx) * t2State(i)
            w(t2idx(i),idx) = w(t2idx(i),idx) + (newState-dw(t2idx(i),idx))*(t-i)
            dw(t2idx(i),idx) = newState
          end do
          idx2t(idx) = t 	! at what time idx was changed
          t2idx(t) = idx 	! what was that idx
        end if
		
        activation = 0
        do j=1,n
          activation = activation + w(j,idx) * state(j)
        end do
      else
        ! Since we only update the upper triangle of w below, so the sum in:
        ! activation = sum(w(:,idx)*state) 
        ! needs to be split to:
        ! activation = sum(w(1:idx,idx)*state(1:idx)) + sum(w(idx,idx+1:n)*state(idx+1:n))
		activation = 0
        do j=1,idx
          activation = activation + w(j,idx)*state(j)
        end do
        do j = idx+1,n
          activation = activation + w(idx,j)*state(j)
        end do
      end if
	  
      if (activation>0) then
        state(idx) = 1
      else
        state(idx) = -1
      end if
	  
      if (dolearn==1) then
        if (speed==1) then
          t2State(t) = state(idx) ! save the state that was changed at time t
        end if
        if (t==1) then
          do j=1,n
            dw(:,j) = state(:)*state(j)
          end do
        else
          if (state(idx)>0) then
            dw(:,idx) = state(:)
            if (speed==0) then
              dw(idx,:) = state(:)
            end if
          else
            dw(:,idx)=-state(:)
            if (speed==0) then
              dw(idx,:) = -state(:)
            end if
          end if
        end if

        if (speed==0) then
          do j=1,n
            w(1:j,j) = w(1:j,j) + dw(1:j,j)
          end do
        end if
        
      end if
      if (t==1) then
        energies(t) = 0
        do j=1,n
          energies(t) = energies(t) - worig(j,j)&
               - 2*state(j)*sum(state(1:j-1)*worig(1:j-1,j))
        end do
      else
        energies(t) = energies(t-1) - 2*(state(idx)-oldState)&
             * (sum(state(1:idx-1) * worig(1:idx-1,idx))&
             + sum(state(idx+1:n) * worig(idx+1:n,idx)))
      end if
    end do
	
    if ((speed==1) .and. (dolearn==1)) then
      t = steps+1
      do idx=1,n
        w(:,idx) = w(:,idx)+dw(:,idx)*(t-idx2t(idx))
        do i = idx2t(idx)+1,t-1
          newState = state(idx)*t2State(i)
          w(t2idx(i),idx) = w(t2idx(i),idx) + (newState-dw(t2idx(i),idx))*(t-i)
          dw(t2idx(i),idx) = newState
        end do
      end do
    end if
  end SUBROUTINE learn
end module sofort
