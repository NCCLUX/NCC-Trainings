#### <u>Shared variable</u>

 - All the threads have access to the shared variable.
 - By default, in the parallel region, all the variables are considered shared variables except the loop iteration counter variables.
 
!!! Note

	Shared variables should be handled carefully; otherwise, they cause race conditions in the program.

<figure markdown>
![](../figures/shared-white.png){align=center width=500}
</figure>

??? example "Examples: Shared variable"

    === "(C/C++)"
    
        ```c
        #include <iostream>
        #include <omp.h>
        
        using namespace std;
        
        int main()
        {
          // Array size
          int N = 10;
          
          // Initialize the variables
          float *a;
          
          // Allocate the memory
          a  = (float*)malloc(sizeof(float) * N);
          
          //#pragma omp parallel for
          // or 
        #pragma omp parallel for shared(a)
          for (int i = 0; i < N; i++)
            {
              a[i] = a[i] + i;  
              cout << "value of a in the parallel region" << a[i] << endl;
            }
            
          for (int i = 0; i < N; i++)
            cout << "value of a after the parallel region " << a[i] << endl;
            
          return 0;
        }
        ```
	
    === "(FORTRAN)"
    
        ```c
        program main
          use omp_lib
          implicit none
          
          ! Input vectors
          real(8), dimension(:), allocatable :: a
          
          integer :: n, i
          n=10
          
          ! Allocate memory for vector
          allocate(a(n))
          
          !$omp parallel shared(a)
          !$omp do
          do i = 1, n
              a(i) = a(i) + i
              print*, 'value of a in the parallel region', a(i)
          end do
          !$omp end do
          !$omp end parallel

          do i = 1, n
              a(i) = a(i) + i
              print*,'value of an after the parallel region', a(i)
          end do

          ! Delete the memory
          deallocate(a)
          
        end program main
        ```

??? question "Questions"

     - Does the value of vector `a` change after the parallel loop, if not why, think?
     - Do we really need to mention `shared(a)`, is it neccessary? 

#### <u>Private variable</u>

 - Each thread will have its own copy of the private variable.
 - And the private variable is only accessible within the parallel region,
 not outside of the parallel region.
 - By default, the loop iteration counters are considered as a private.
 - A change made by one thread is not visible to other threads.


<figure markdown>
![](../figures/private-white.png){align=center width=500}
</figure>


??? example "Examples: Private variable"

    === "(C/C++)"
    
        ```c
        #include <iostream>
        #include <omp.h>
        
        using namespace std;
        
        int main()
        {
          // Array size
          int N = 10;
          
          // Initialize the variables
          float *a,b,c;
          b = 1.0;
          c = 2.0;
	  
          // Allocate the memory
          a  = (float*)malloc(sizeof(float) * N);
          
        #pragma omp parallel for private(b,c)
          for (int i = 0; i < N; i++)
            {
              b = a[i] + i;
              c = b + 10 * i;
              cout << "value of c in the parallel region " << c << endl;
            }
          
          cout << "value of c after the parallel region " << c << endl;	
                        
          return 0;
        }
        ```	

    === "(FORTRAN)"
    
        ```c
        program main
          use omp_lib
          implicit none
          
          ! Input vectors
          real(8), dimension(:), allocatable :: a
          
          real(8) :: b, c
          integer :: n, i  
          n=10
          b=1.0
          c=2.0
          
          ! Allocate memory for vector
          allocate(a(n))
          
          !$omp parallel private(b,c) shared(a)
          !$omp do
          do i = 1, n
              b = a(i) + i
              c = b + 10 * i
              print*, 'value of c in the parallel region', c
          end do
          !$omp end do
          !$omp end parallel
          print*, 'value of c after the parallel region', c
	  
          ! Delete the memory
          deallocate(a)
          
        end program main
        ```

??? question "Questions"

     - What is the value of the variable `c` in the parallel region and after the parallel region?
     - After the parallel region, has variable `c` been updated or not? 


#### <u>Lastprivate</u>

 - lastprivate: is also similar to a private clause
 - But each thread will have an uninitialized copy of the variables passed
 as lastprivate
 - At the end of the parallel loop or sections, the final variable value will
 be the last thread accessed value in the section or in a parallel loop.

??? example "Examples: Lastprivate variable"

    === "(C/C++)"
    
        ```c
        #include<iostream>
        #include<omp.h>
        
        using namespace std;

        int main()
        {
          int n = 10;
          int var = 5;
          omp_set_num_threads(10);
        #pragma omp parallel for lastprivate(var)
          for(int i = 0; i < n; i++)
            {
              var += omp_get_thread_num();
              cout << " lastprivate in the parallel region " << var << endl;
            } /*-- End of parallel region --*/
          cout << "lastprivate after the parallel region " << var <<endl;
  
          return 0;
        }
        ```	


    === "(FORTRAN)"
    
        ```c
        program main
          use omp_lib
          implicit none
  
          ! Initialise the variable
          real(8) :: var
          integer :: n, i  
          n = 10
          var = 5
  
          call omp_set_num_threads(10)
  
          !$omp parallel 
          !$omp do lastprivate(var)
          do i = 1, n
             var  =  var + omp_get_thread_num()
             print*, 'lastprivate in the parallel region ', var
          end do
          !$omp end do
          !$omp end parallel

          print*, 'lastprivate after the parallel region ', var
 
        end program main
        ```

??? question "Questions"

     - What is the value of the variable `var` in the parallel region and after the parallel region?
     - Do you think the initial value of variable `var` is considered within the parallel region? 


#### <u>Firstprivate</u>

 - firstprivate: is similar to a private clause
 - But each thread will have an initialized copy of the variables passed
 as firstprivate
 - Available for parallel constructs, loops, sections and single
 constructs

??? example "Examples: Firstprivate variable"

    === "(C/C++)"
    
        ```c
        #include<iostream>
        #include<omp.h>
        
        using namespace std;

        int main()
        {
          int n = 10;
          int var = 5;
          omp_set_num_threads(10);
        #pragma omp parallel for firstprivate(var)
          for(int i = 0; i < n; i++)
            {
              var += omp_get_thread_num();
              cout << " lastprivate in the parallel region " << var << endl;
            } /*-- End of parallel region --*/
          cout << "lastprivate after the parallel region " << var <<endl;
  
          return 0;
        }
        ```	


    === "(FORTRAN)"
    
        ```c
        program main
          use omp_lib
          implicit none
  
          ! Initialise the variable
          real(8) :: var
          integer :: n, i  
          n = 10
          var = 5
  
          call omp_set_num_threads(10)
  
          !$omp parallel 
          !$omp do firstprivate(var)
          do i = 1, n
             var  =  var + omp_get_thread_num()
             print*, 'lastprivate in the parallel region ', var
          end do
          !$omp end do
          !$omp end parallel

          print*, 'lastprivate after the parallel region ', var
 
        end program main
        ```


??? question "Questions"

     - What is the value of the variable `var` in the parallel region and after the parallel region?
     - Is variable `var` has been updated after the parallel region, if not why, think?


