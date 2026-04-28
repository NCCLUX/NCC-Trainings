#### <u>Parallel Construct</u>

In this exercise, we will create a parallel region and execute the computational content in parallel.
First, however, this exercise is to create a parallel region and understand the threads' behaviour in parallel.
In later exercises, we will study how to parallelise the computational task within the parallel region.

<figure markdown>
![](../figures/parallel-white.png){align=center width=500}
</figure>

To create a parallel region, we use the following parallel constructs:

!!! Info "Parallel Constructs"

    === "C/C++"
        ```
        #pragma omp parallel
        ```
    === "FORTRAN"
        ```
        !$omp parallel 
        ```

The above figure illustrates the parallel region behaviour;
as we notice, within the parallel region, we get parallel threads.
This means parallel threads can be executed independently of each other,
and there is no order of execution.

At the same time, in order to enable OpenMP constructs, clauses, and environment variables. etc., we need to include the OpenMP library as follows:

!!! Info "OpenMP library"

    === "C/C++"
        ```
        #include<omp.h>
        ```
    === "FORTRAN"
        ```
        use omp_lib
        ```


#### <u>Compilers</u>

The following compilers would support the OpenMP programming model.

 - [GNU](https://gcc.gnu.org/) - It is an open source and can be used for Intel and AMD CPUs
 - [Intel](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html#gs.zd201n) - It is from Intel and only optimized for Intel CPUs
 - [AOOC](https://www.amd.com/content/dam/amd/en/documents/pdfs/developer/aocc/aocc-v4.0-ga-user-guide.pdf) - Suitable for AMD CPUs, especially “Zen” core architecture.


!!! Info "Examples (GNU, Intel and AMD): Compilation"

    === "GNU"
        ```c
        $ gcc test.c -fopenmp
        $ g++ test.cc -fopenmp
        $ gfortran test.f90 -fopenmp
        ```
    === "Intel"
        ```c
        $ icc test.c -qopenmp
        $ icpc test.cc -qopenmp
        $ ifort test.f90 -qopenmp        
        ```
    === "AOOC"
        ```c
        $ clang test.c -fopenmp
        $ clang++ test.cc -fopenmp
        $ flang test.f90 -fopenmp        
        ```

### <u>Questions and Solutions</u>


??? Example "Examples: Hello World"

    === "Serial-version (C/C++)"
        ``` c
        #include<iostream>
        using namespace std;
        
        int main()
        {
          cout << endl;
          cout << "Hello world from the master thread"<< endl;
          cout << endl;
          
          return 0;
        }
        ```

    === "Serial-version (FORTRAN)"
        ``` fortran
        program Hello_world_Serial
        
        print *, 'Hello world from master thread'
        
        end program
        ```
	
    === "OpenMP-version (C/C++)"
        ``` c
        #include<iostream>
        #include<omp.h>
        using namespace std;
        int main()
        {
          cout << "Hello world from the master thread "<< endl;
          cout << endl;
          
          // creating the parallel region (with N number of threads)
          #pragma omp parallel
           {
                cout << "Hello world from the parallel region "<< endl;
            } // parallel region is closed
            
        cout << endl;
        cout << "end of the programme from the master thread" << endl;
        return 0;
        }
        ```

    === "OpenMP-version (FORTRAN)"
        ``` fortran
        program Hello_world_OpenMP
        use omp_lib

        print *, 'Hello world from master thread'

        !$omp parallel
        print *, 'Hello world from parallel region'
        !$omp end parallel

        print *,'end of the programme from master thread'
        
        end program
        ```

??? "Compilation and Output"

    === "Serial-version (C/C++)"
        ```c
        // compilation
        $ g++ Hello-world-Serial.cc -o Hello-World-Serial
        
        // execution 
        $ ./Hello-World-Serial
        
        // output
        $ Hello world from the master thread
        ```

    === "Serial-version (FORTRAN)"
        ```c
        // compilation
        $ gfortran Hello-world-Serial.f90 -o Hello-World-Serial
        
        // execution 
        $ ./Hello-World-Serial
        
        // output
        $ Hello world from the master thread
        ```
        
    === "OpenMP-version (C/C++)"
        ```c
        // compilation
        $ g++ -fopenmp Hello-world-OpenMP.cc -o Hello-World-OpenMP
        
        // execution
        $ ./Hello-World-OpenMP
        
        // output
        $ Hello world from the parallel region
        Hello world from the parallel region
        ..
        ..
        Hello world from the parallel region
        
        end of the programme from the master thread
        ```

    === "OpenMP-version (FORTRAN)"
        ```c
        // compilation
        $ gfortran -fopenmp Hello-world-OpenMP.f90 -o Hello-World-OpenMP
        
        // execution
        $ ./Hello-World-OpenMP
        
        // output
        $ Hello world from the master thread
        Hello world from the parallel region
        ..
        ..
        Hello world from the parallel region
        End of the programme from the master thread
        ```

??? question "Questions"

     - What do you notice from those examples? Can you control parallel region printout, that is, how many times it should be printed or executed?     
     - What happens if you do not use the OpenMP library, `#include<omp.h> or use omp_lib`?


Although creating a parallel region would allow us to do the parallel computation, however, at the same time, we should have control over the threads being created in the parallel region, for example, how many threads are needed for a particular computation, thread number, etc. For this, we need to know a few of the important environmental routines provided by OpenMP. The list below shows a few of the most important environment routines that the programmer should know about for optimised OpenMP coding.

#### <u>Environment Routines (important)</u>

 - Define the number of threads to be used within the parallel region
 
        (C/C++): void omp_set_num_threads(int num_threads);
        (FORTRAN): subroutine omp_set_num_threads(num_threads) 
        integer num_threads

 - To get the number of threads in the current parallel region

        (C/C++): int omp_get_num_threads(void);
        (FORTRAN): integer function omp_get_num_threads()

 - To get available maximum threads (system default)

        (c/c++): int omp_get_max_threads(void);
        (FORTRAN): integer function omp_get_max_threads()

 - To get thread numbers (e.g., 1, 4, etc.)

        (c/c+): int omp_get_thread_num(void);
        (FORTRAN): integer function omp_get_thread_num()

 - To know the number of processors available to the device
 
        (c/c++): int omp_get_num_procs(void);
        (FROTRAN): integer function omp_get_num_procs()


### <u>Questions and Solutions</u>


??? question "Questions"

     - How can you identify the thread numbers within the parallel region?
     - What happens if you not set `omp_set_num_threads()`, for example, `omp_set_num_threads(5)|call omp_set_num_threads(5)`, what do you notice? 
     Alternatively, you can also set the number of threads to be used in the application during the compilation `export OMP_NUM_THREADS`; what do you see?

    === "Question (C/C++)"

        ```c
        #include<iostream>
        #include<omp.h>
        using namespace std;
        int main()
        {
          cout << "Hello world from the master thread "<< endl;
          cout << endl;
                    
          // creating the parallel region (with N number of threads)
          #pragma omp parallel
           {
                //cout << "Hello world from thread id "
                << " from the team size of "
                << endl;
            } // parallel region is closed
            
        cout << endl;
        cout << "end of the programme from the master thread" << endl;
        return 0;
        }
        ```

    === "Question (FORTRAN)"
        ``` fortran
        program Hello_world_OpenMP
        use omp_lib
                
        !$omp parallel 
        !! print *, 
        !$omp end parallel
        
        end program
        ```
	
    === "Answer (C/C++)"
        ```c
        #include<iostream>
        #include<omp.h>
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

    === "Answer (FORTRAN)"
        ``` fortran
        program Hello_world_OpenMP
        use omp_lib
                
        !$omp parallel 
        print *, 'Hello world from thread id ', omp_get_thread_num(), 'from the team size of', omp_get_num_threads()
        !$omp end parallel
        
        end program
        ```
        
    === "Answer"
        ``` c
        $ export OMP_NUM_THREADS=10
        // or 
        $ setenv OMP_NUM_THREADS 4
        // or
        $ OMP NUM THREADS=4 ./omp code.exe
        ```

    === "Solution Output (C/C++)"
        ```c
        ead id Hello world from thread id Hello world from thread id 3 from the team size of 9 from the team size of 52 from the team size of  from the team size of 10
        0 from the team size of 10
        10
        10
        10
        7 from the team size of 10
        4 from the team size of 10
        8 from the team size of 10
        1 from the team size of 10
        6 from the team size of 10
        ```
        
    === "Solution Output (FORTRAN)"

        ```c
        Hello world from thread id            0 from the team size of          10
        Hello world from thread id            4 from the team size of          10
        Hello world from thread id            5 from the team size of          10
        Hello world from thread id            9 from the team size of          10
        Hello world from thread id            2 from the team size of          10
        Hello world from thread id            3 from the team size of          10
        Hello world from thread id            7 from the team size of          10
        Hello world from thread id            6 from the team size of          10
        Hello world from thread id            8 from the team size of          10
        Hello world from thread id            1 from the team size of          10
        ```


#### <u>Utilities</u>


The main aim is to do the parallel computation to speed up computation on a given parallel architecture. Therefore, measuring the timing and comparing the solution between serial and parallel code is very important. In order to measure the timing, OpenMP provides an environmental variable, `omp_get_wtime()`.

??? Info "Time measuring"

    === "C/C++"
        ```
        double start; 
        double end; 
        start = omp_get_wtime(); 
        ... work to be timed ... 
        end = omp_get_wtime(); 
        printf("Work took %f seconds\n", end - start);
        ```
	
    === "FORTRAN"
        ```
        DOUBLE PRECISION START, END 
        START = omp_get_wtime() 
        ... work to be timed ... 
        END = omp_get_wtime() 
        PRINT *, "Work took", END - START, "seconds"        
        ```


