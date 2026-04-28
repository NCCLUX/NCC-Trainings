### <u>[Loop scheduling](https://www.openmp.org/spec-html/5.0/openmpsu41.html#x64-1370002.9.2.1)</u>

   However, the above example is very simple.
   Because, in most cases, we would end up doing a large list of arrays with complex computations within the loop.
   Therefore, the workloading should be optimally distributed among the threads in those cases.
   To handle those considerations, OpenMP has provided the following loop-sharing clauses. They are: `Static`, `Dynamic`, `Guided`, `Auto`, and `Runtime`.

<figure markdown>
![](../figures/schedule.png){width="500", align=middle}
<figcaption></figcaption>
</figure>

??? Info "Example - Loop scheduling clauses"

    === "C/C++"

        ```c
        #pragma omp parallel for schedule(static)
        for(int i = 0; i < n; i ++)
          {
            c[i] = a[i] + b[i];
          }
          
        //or
        
        #pragma omp parallel 
        #pragma omp for schedule(static)
        for(int i = 0; i < n; i ++)
          {
            c[i] = a[i] + b[i];
          }
        ```  
        
    === "FORTRAN"
    
        ```c
        !$omp parallel do schedule(static)
        do i = 1, n
          c(i) = a(i) + b(i)
        end do
        !$omp end parallel do
        
        //or
        
        !$omp parallel
        !$omp do schedule(static)
        do i = 1, n
          c(i) = a(i) + b(i)
        end do
        !$omp end do
        !$omp end parallel
        ```
	
#### <u>Static</u>

 - The number of iterations are divided by chunksize. 
 - If the chunksize is not provided, a number of iterations will be divided by the size of the team of threads.
    - e.g., n=100, numthreads=5; each thread will execute the 20 iterations in parallel.
 - This is useful when the computational cost is similar to each iteration.

??? example "Examples and Question: static"

    === "OpenMP(C/C++)"
        ```c
        #include <iostream>
        #include <omp.h>
        using namespace std;
        int main()
        {
         int N = 10;
         omp_set_num_threads(5);
         
        #pragma omp parallel for schedule(static)
        for(int i = 0; i < N; i++)
           {
            cout << " Thread id" << " " << omp_get_thread_num() << endl;    
           }  
          return 0;
        }
        ```


    === "OpenMP(FORTRAN)"
        
        ```c
        program main
        use omp_lib
        implicit none
        
        integer :: n, i  
        n = 10
        
        call omp_set_num_threads(5)

        !$omp parallel
        !$omp do schedule(static)
        do i = 1, n
          print*, 'Thread id', omp_get_thread_num()
        end do
        !$omp end do
        !$omp end parallel
        
        end program main
        ```
        
    === "Output"
        ```c
        Thread id           0
        Thread id           0
        Thread id           4
        Thread id           4
        Thread id           3
        Thread id           3
        Thread id           2
        Thread id           2
        Thread id           1
        Thread id           1
        ```

    - What happens if you would set the chunksize, for example, `schedule(static,4)`? What do you notice?
 
#### <u>Dynamic</u>

 - The number of iterations are divided by chunksize.
 - If the chunksize is not provided, the default value will be considered 1.
 - This is useful when the computational cost is different in the iteration.
 - This will quickly place the chunk of data in the queue.

??? example "Examples and Question: dynamic"

    === "OpenMP(C/C++)"
        ```c
        #include <iostream>
        #include <omp.h>
        using namespace std;
        int main()
        {
         int N = 10;
         omp_set_num_threads(5);
         
        #pragma omp parallel for schedule(dynamic)
        for(int i = 0; i < N; i++)
           {
            cout << " Thread id" << " " << omp_get_thread_num() << endl;    
           }  
          return 0;
        }
        ```


    === "OpenMP(FORTRAN)"
        
        ```c
        program main
        use omp_lib
        implicit none
        
        integer :: n, i  
        n = 10
        
        call omp_set_num_threads(5)

        !$omp parallel
        !$omp do schedule(dynamic)
        do i = 1, n
          print*, 'Thread id', omp_get_thread_num()
        end do
        !$omp end do
        !$omp end parallel
        
        end program main
        ```
        
    === "Output"
        ```c
        Thread id  Thread id 20 Thread id 
        4 Thread id 2
        Thread id 2
        Thread id 2
        Thread id 2
        Thread id 2
        Thread id
        Thread id 1 
        3
        ```
    - What happens if you would set the chunksize, for example, schedule(dynamic,4)? What do you notice?
    - Do you notice if the iterations are divided by the chunksize that we set?
    
#### <u>Guided</u>

 - Similar to dynamic scheduling, the number of iterations are divided by chunksize.
 - But the chunk of the data size is decreasing, which is proportional to the number of unsigned iterations divided by the number of threads.
 - If the chunksize is not provided, the default value will be considered 1.
 - This is useful when there is poor load balancing at the end of the iteration.

??? example "Examples and Question: guided"

    === "OpenMP(C/C++)"
        ```c
        #include <iostream>
        #include <omp.h>
        using namespace std;
        int main()
        {
         int N = 10;
         omp_set_num_threads(5);
         
        #pragma omp parallel for schedule(guided)
        for(int i = 0; i < N; i++)
           {
            cout << " Thread id" << " " << omp_get_thread_num() << endl;    
           }  
          return 0;
        }
        ```


    === "OpenMP(FORTRAN)"
        
        ```c
        program main
        use omp_lib
        implicit none
        
        integer :: n, i  
        n = 10
        
        call omp_set_num_threads(5)

        !$omp parallel
        !$omp do schedule(guided)
        do i = 1, n
          print*, 'Thread id', omp_get_thread_num()
        end do
        !$omp end do
        !$omp end parallel
        
        end program main
        ```
        
    === "Output"
        ```c
        Thread id Thread id   Thread id0 41
        Thread id 0
        
        Thread id 4
        Thread id 4
        Thread id 2
        Thread id 2
        Thread id 3 Thread id
        ```

    - Are there any differences between `auto` and `guided` or `dynamic`?

#### <u>Auto</u>

 - Here the compiler chooses the best combination of the chunksize to be used. 

??? example "Examples and Question: auto"

    === "OpenMP(C/C++)"
        ```c
        #include <iostream>
        #include <omp.h>
        using namespace std;
        int main()
        {
         int N = 10;
         omp_set_num_threads(5);
         
        #pragma omp parallel for schedule(auto)
        for(int i = 0; i < N; i++)
           {
            cout << " Thread id" << " " << omp_get_thread_num() << endl;    
           }  
          return 0;
        }
        ```


    === "OpenMP(FORTRAN)"
        
        ```c
        program main
        use omp_lib
        implicit none
        
        integer :: n, i  
        n = 10
        
        call omp_set_num_threads(5)

        !$omp parallel
        !$omp do schedule(auto)
        do i = 1, n
          print*, 'Thread id', omp_get_thread_num()
        end do
        !$omp end do
        !$omp end parallel
        
        end program main
        ```
        
    === "Output"
        ```c
        Thread id Thread id Thread id    Thread id0 34 Thread id 
        Thread id 0
        1
         Thread id 1

         Thread id 3
        2
         Thread id 2
        
         Thread id 4
        ```
    - What would you choose for your application, auto, dynamic, guided, or static? If you are going to choose either one of them, then have a valid reason. 


#### <u>[Runtime](https://www.openmp.org/spec-html/5.0/openmpse49.html#x288-20520006.1)</u>

  - During the compilation, we simply set the loop scheduling concept.


??? Info "Example: Loop scheduling clauses - runtime"

    === "Compilation"

        ```c
        setenv OMP_SCHEDULE="guided,4" 
        setenv OMP_SCHEDULE="dynamic" 
        setenv OMP_SCHEDULE="nonmonotonic:dynamic,4"
        // or
        export OMP_SCHEDULE="guided,4" 
        export OMP_SCHEDULE="dynamic" 
        export OMP_SCHEDULE="nonmonotonic:dynamic,4"
        ```
	
??? example "Examples and Question: runtime"

    === "OpenMP(C/C++)"
        ```c
        #include <iostream>
        #include <omp.h>
        using namespace std;
        int main()
        {
         int N = 10;
         omp_set_num_threads(5);
         
        #pragma omp parallel for schedule(runtime)
        for(int i = 0; i < N; i++)
           {
            cout << " Thread id" << " " << omp_get_thread_num() << endl;    
           }  
          return 0;
        }
        ```

    === "OpenMP(FORTRAN)"
        
        ```c
        program main
        use omp_lib
        implicit none
        
        integer :: n, i  
        n = 10
        
        call omp_set_num_threads(5)

        !$omp parallel
        !$omp do schedule(runtime)
        do i = 1, n
          print*, 'Thread id', omp_get_thread_num()
        end do
        !$omp end do
        !$omp end parallel
        
        end program main
        ```
	
    === "Compilation"
        
        ```c
        export OMP_SCHEDULE="dynamic,3"
        // check if you have exported the environment value
        $ env | grep OMP_SCHEDULE
        $ OMP_SCHEDULE=dynamic,3 
        // if you want to unset
        $ unset OMP_SCHEDULE
        $ env | grep OMP_SCHEDULE
        // it(OMP_SCHEDULE=dynamic,3) will be removed
        ```
