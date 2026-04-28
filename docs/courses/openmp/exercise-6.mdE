<figure markdown>
![](../figures/simd.png){align=center width=500}
<figcaption></figcaption>
</figure>


In this exercise, we will try to add the `simd` classes to our existing problems, for example, vector addition. 


??? example "Examples and Question: SIMD - Vector Addition"


    === "Serial(C/C++)"
    
        ```c  
        #include <stdio.h>
        #include <stdlib.h>
        #include <math.h>
        #include <assert.h>
        #include <time.h>
        
        #define N 5120
        #define MAX_ERR 1e-6

        // CPU function that adds two vector 
        float * Vector_Add(float *a, float *b, float *c, int n) 
        {
          for(int i = 0; i < n; i ++)
            {
              c[i] = a[i] + b[i];
            }
          return c;
        }

        int main()
        {
          // Initialize the variables
          float *a, *b, *c;       
  
          // Allocate the memory
          a   = (float*)malloc(sizeof(float) * N);
          b   = (float*)malloc(sizeof(float) * N);
          c = (float*)malloc(sizeof(float) * N);
  
          // Initialize the arrays
          for(int i = 0; i < N; i++)
            {
              a[i] = 1.0f;
              b[i] = 2.0f;
            }
    
          // Start measuring time
          clock_t start = clock();

          // Executing vector addition function 
          Vector_Add(a, b, c, N);

          // Stop measuring time and calculate the elapsed time
          clock_t end = clock();
          double elapsed = (double)(end - start)/CLOCKS_PER_SEC;
        
          printf("Time measured: %.3f seconds.\n", elapsed);
  
          // Verification
          for(int i = 0; i < N; i++)
            {
              assert(fabs(c[i] - a[i] - b[i]) < MAX_ERR);
            }

          printf("c[0] = %f\n", c[0]);
          printf("PASSED\n");
    
          // Deallocate the memory
          free(a); 
          free(b); 
          free(c);
   
          return 0;
        }
        ```

    === "Serial(FORTRAN)"
        ```c
        module Vector_Addition_Mod  
        implicit none 
          contains
        subroutine Vector_Addition(a, b, c, n)
        ! Input vectors
        real(8), intent(in), dimension(:) :: a
        real(8), intent(in), dimension(:) :: b
        real(8), intent(out), dimension(:) :: c
        integer :: i, n
          do i = 1, n
            c(i) = a(i) + b(i)
          end do
         end subroutine Vector_Addition
        end module Vector_Addition_Mod

        program main
        use Vector_Addition_Mod
        implicit none
        ! Input vectors
        real(8), dimension(:), allocatable :: a
        real(8), dimension(:), allocatable :: b 
        ! Output vector
        real(8), dimension(:), allocatable :: c
        ! real(8) :: sum = 0

        integer :: n, i  
        print *, "This program does the addition of two vectors "
        print *, "Please specify the vector size = "
        read *, n

        ! Allocate memory for vector
        allocate(a(n))
        allocate(b(n))
        allocate(c(n))
  
        ! Initialize content of input vectors, 
        ! vector a[i] = sin(i)^2 vector b[i] = cos(i)^2
        do i = 1, n
          a(i) = sin(i*1D0) * sin(i*1D0)
          b(i) = cos(i*1D0) * cos(i*1D0) 
        enddo
    
        ! Call the vector addition subroutine 
        call Vector_Addition(a, b, c, n)

        !!Verification
        do i = 1, n
          if (abs(c(i)-(a(i)+b(i)) == 0.00000)) then 
           else
             print *, "FAIL"
           endif
        enddo
        print *, "PASS"
    
        ! Delete the memory
        deallocate(a)
        deallocate(b)
        deallocate(c)
  
        end program main

        ```



    === "Template(C/C++)"
    
        ```c
        #include <stdio.h>
        #include <stdlib.h>
        #include <math.h>
        #include <assert.h>
        #include <time.h>
        #include <omp.h>
        
        #define N 5120
        #define MAX_ERR 1e-6

        // CPU function that adds two vector 
        float * Vector_Add(float *a, float *b, float *c, int n) 
        {
        // ADD YOUR PARALLEL REGION FOR THE LOOP SIMD
          for(int i = 0; i < n; i ++)
            {
              c[i] = a[i] + b[i];
            }
          return c;
        }

        int main()
        {
          // Initialize the variables
          float *a, *b, *c;       
  
          // Allocate the memory
          a   = (float*)malloc(sizeof(float) * N);
          b   = (float*)malloc(sizeof(float) * N);
          c = (float*)malloc(sizeof(float) * N);
  
          // Initialize the arrays
          for(int i = 0; i < N; i++)
            {
              a[i] = 1.0f;
              b[i] = 2.0f;
            }
    
          // Start measuring time
          clock_t start = clock();

          // ADD YOUR PARALLEL REGION HERE	
          // Executing vector addition function 
          Vector_Add(a, b, c, N);

          // Stop measuring time and calculate the elapsed time
          clock_t end = clock();
          double elapsed = (double)(end - start)/CLOCKS_PER_SEC;
        
          printf("Time measured: %.3f seconds.\n", elapsed);
  
          // Verification
          for(int i = 0; i < N; i++)
            {
              assert(fabs(c[i] - a[i] - b[i]) < MAX_ERR);
            }

          printf("c[0] = %f\n", c[0]);
          printf("PASSED\n");
    
          // Deallocate the memory
          free(a); 
          free(b); 
          free(c);
   
          return 0;
        }

        ```
        
    === "Template(FORTRAN)"
        ```c
        module Vector_Addition_Mod  
        implicit none 
          contains
        subroutine Vector_Addition(a, b, c, n)
        use omp_lib
        ! Input vectors
        real(8), intent(in), dimension(:) :: a
        real(8), intent(in), dimension(:) :: b
        real(8), intent(out), dimension(:) :: c
        integer :: i, n
        !! ADD YOUR PARALLEL DO LOOP WITH SIMD
          do i = 1, n
            c(i) = a(i) + b(i)
          end do
         end subroutine Vector_Addition
        end module Vector_Addition_Mod

        program main
        use Vector_Addition_Mod
        implicit none
        ! Input vectors
        real(8), dimension(:), allocatable :: a
        real(8), dimension(:), allocatable :: b 
        ! Output vector
        real(8), dimension(:), allocatable :: c
        ! real(8) :: sum = 0

        integer :: n, i  
        print *, "This program does the addition of two vectors "
        print *, "Please specify the vector size = "
        read *, n

        ! Allocate memory for vector
        allocate(a(n))
        allocate(b(n))
        allocate(c(n))
  
        ! Initialize content of input vectors, 
        ! vector a[i] = sin(i)^2 vector b[i] = cos(i)^2
        do i = 1, n
          a(i) = sin(i*1D0) * sin(i*1D0)
          b(i) = cos(i*1D0) * cos(i*1D0) 
        enddo

        !! ADD YOUR PARALLEL REGION 
        ! Call the vector add subroutine 
        call Vector_Addition(a, b, c, n)

        !!Verification
        do i = 1, n
          if (abs(c(i)-(a(i)+b(i)) == 0.00000)) then 
           else
             print *, "FAIL"
           endif
        enddo
        print *, "PASS"
    
        ! Delete the memory
        deallocate(a)
        deallocate(b)
        deallocate(c)
  
        end program main

        ```

    === "Solution(C/C++)"
    
        ```c       
        #include <stdio.h>
        #include <stdlib.h>
        #include <math.h>
        #include <assert.h>
        #include <time.h>
        #include <omp.h>
        
        #define N 5120
        #define MAX_ERR 1e-6

        // CPU function that adds two vector 
        float * Vector_Add(float *a, float *b, float *c, int n) 
        {
        // ADD YOUR PARALLEL SIMD
        #pragma omp for simd
          for(int i = 0; i < n; i ++)
            {
              c[i] = a[i] + b[i];
            }
          return c;
        }

        int main()
        {
          // Initialize the variables
          float *a, *b, *c;       
  
          // Allocate the memory
          a   = (float*)malloc(sizeof(float) * N);
          b   = (float*)malloc(sizeof(float) * N);
          c = (float*)malloc(sizeof(float) * N);
  
          // Initialize the arrays
          for(int i = 0; i < N; i++)
            {
              a[i] = 1.0f;
              b[i] = 2.0f;
            }
    
          double start = omp_get_wtime();
          #pragma omp parallel 
          // Executing vector addition function 
          Vector_Add(a, b, c, N);
          double end = omp_get_wtime();
          printf("Work took %f seconds\n", end - start);
  
          // Verification
          for(int i = 0; i < N; i++)
            {
              assert(fabs(c[i] - a[i] - b[i]) < MAX_ERR);
            }

          printf("c[0] = %f\n", c[0]);
          printf("PASSED\n");
    
          // Deallocate the memory
          free(a); 
          free(b); 
          free(c);
   
          return 0;
        }

        ```

    === "Solution(FORTRAN)"
        ```c
        module Vector_Addition_Mod  
        implicit none 
          contains
        subroutine Vector_Addition(a, b, c, n)
        use omp_lib
        ! Input vectors
        real(8), intent(in), dimension(:) :: a
        real(8), intent(in), dimension(:) :: b
        real(8), intent(out), dimension(:) :: c
        integer :: i, n
        !$omp do simd
          do i = 1, n
            c(i) = a(i) + b(i)
          end do
        !$omp end do simd
         end subroutine Vector_Addition
        end module Vector_Addition_Mod

        program main
        use Vector_Addition_Mod
        implicit none
        ! Input vectors
        real(8), dimension(:), allocatable :: a
        real(8), dimension(:), allocatable :: b 
        ! Output vector
        real(8), dimension(:), allocatable :: c
        ! real(8) :: sum = 0
        double precision :: start, end
	
        integer :: n, i  
        print *, "This program does the addition of two vectors "
        print *, "Please specify the vector size = "
        read *, n

        ! Allocate memory for vector
        allocate(a(n))
        allocate(b(n))
        allocate(c(n))
  
        ! Initialize content of input vectors, 
        ! vector a[i] = sin(i)^2 vector b[i] = cos(i)^2
        do i = 1, n
          a(i) = sin(i*1D0) * sin(i*1D0)
          b(i) = cos(i*1D0) * cos(i*1D0) 
        enddo

        start = omp_get_wtime()
        !$omp parallel 
        ! Call the vector addition subroutine 
        call Vector_Addition(a, b, c, n)
        !$omp end parallel
        end = omp_get_wtime()
        PRINT *, "Work took", end - start, "seconds"
	
        !!Verification
        do i = 1, n
          if (abs(c(i)-(a(i)+b(i)) == 0.00000)) then 
           else
             print *, "FAIL"
           endif
        enddo
        print *, "PASS"
        
        ! Delete the memory
        deallocate(a)
        deallocate(b)
        deallocate(c)
        
        end program main
        ```
       - Please try the examples without the `simd` clause. Do you notice any performance differences? 




####<u>Critical, Single, and Master</u>

We will explore how single, master and critical are working in the OpenMP programming model. For this, we consider the following simple examples.


??? "Examples and Question: Critical, Single and Master"

    === "(C/C++)"
        ```c
        #include<iostream>
        #include<omp.h>
        #include <omp.h>
        using namespace std;
        int main()
        {
          cout << "Hello world from the master thread "<< endl;
          cout << endl;
                    
          // creating the parallel region (with N number of threads)
          #pragma omp parallel
           {
                cout << "Hello world from thread id "
                << omp_get_thread_num() << " from the team size of "
                << omp_get_num_threads()
                << endl;
            } // parallel region is closed
            
        cout << endl;
        cout << "end of the programme from the master thread" << endl;
        return 0;
        }
        ```

    === "FORTRAN)"
        ``` fortran
        program Hello_world_OpenMP
        use omp_lib
                
        !$omp parallel 
        print *, 'Hello world from thread id ', omp_get_thread_num(), 'from the team size of', omp_get_num_threads()
        !$omp end parallel
        
        end program
        ```


      - Try single clause
      - Try master clause
      - Try critical clause

