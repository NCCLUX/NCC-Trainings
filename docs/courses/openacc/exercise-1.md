#### <u>Compute Constructs</u>


In this initial exercise, we will explore the process of offloading computational tasks to a device, specifically the GPU (Graphics Processing Unit). The primary objective of OpenACC is to streamline this offloading using its dedicated APIs.

OpenACC offers two primary constructs for offloading computations to the GPU, which we will discuss in detail:

1. **Parallel Construct (`parallel`)**: This construct allows for the parallelization of computations across multiple processing units. It is a suitable choice for programmers who possess a strong understanding of their code's parallel behavior.

2. **Kernels Construct (`kernels`)**: This alternative also enables parallelization but provides greater control over the parallel region. It is generally recommended for programmers who may not be as familiar with the intricacies of parallel execution, as the compiler will manage the complexities of safe parallelization within this construct.

Both constructs serve similar purposes in facilitating computation on the GPU; however, the choice between them depends on the programmer's comfort level with parallel programming. If you have significant experience and knowledge of the computations being executed in parallel, you may opt for the `parallel` construct. In contrast, for those who prefer a more managed approach to ensure safety and correctness in parallel execution, using `kernels` is advisable.

To effectively utilize OpenACC constructs, clauses, and environment variables, it is essential to include the OpenACC library in your code. This inclusion enables access to the full range of OpenACC features and functionalities.

!!! Info "OpenACC library"

    === "C/C++"
        ```
        #include<openacc.h>
        ```
    === "FORTRAN"
        ```
        use openacc
        ```

<figure markdown>
![](../figures/acc-parallel.png){align=center width=500}
</figure>


To create a parallel region in OpenACC, we utilize the following compute constructs:

!!! Info "Parallel Constructs"

    === "C/C++"
        ```
        #pragma acc parallel [clause-list] new-line
           structured block
        ```
    === "FORTRAN"
        ```
        !$acc parallel [ clause-list ]
            structured block
        !$acc end parallel
        ```

??? "Available clauses for parallel"

    === "C/C++ and FORTRAN"
	```c
        async [ ( int-expr ) ]
        wait [ ( int-expr-list ) ]
        num_gangs( int-expr )
        num_workers( int-expr )
        vector_length( int-expr )
        device_type( device-type-list )
        if( condition )
        self [ ( condition ) ]
        reduction( operator : var-list )
        copy( var-list )
        copyin( [ readonly: ] var-list )
        copyout( [ zero: ] var-list )
        create( [ zero: ] var-list )
        no_create( var-list )
        present( var-list )
        deviceptr( var-list )
        attach( var-list )
        private( var-list )
        firstprivate( var-list )
        default( none | present )
	```



!!! Info "Kernels Constructs"

    === "C/C++"
        ```
        #pragma acc kernels [ clause-list ] new-line
           structured block
        ```
	
    === "FORTRAN"
        ```
        !$acc kernels [ clause-list ]
           structured block
        !$acc end kernels
        ```


??? "Available clauses for kernels"

    === "C/C++ and FORTRAN"
        ```c
        async [ ( int-expr ) ]
        wait [ ( int-expr-list ) ]
        num_gangs( int-expr )
        num_workers( int-expr )
        vector_length( int-expr )
        device_type( device-type-list )
        if( condition )
        self [ ( condition ) ]
        copy( var-list )
        copyin( [ readonly: ] var-list )
        copyout( [ zero: ] var-list )
        create( [ zero: ] var-list )
        no_create( var-list )
        present( var-list )
        deviceptr( var-list )
        attach( var-list )
        default( none | present )
        ```

#### <u>Compilers Supporting OpenACC Programming Model</u>

The following compilers provide support for the [OpenACC programming model](https://www.openacc.org/tools), which facilitates the development of parallel applications across various architectures:

- **[GNU Compiler Collection (GCC)](https://gcc.gnu.org/)**: This is an open-source compiler that supports both Nvidia and AMD CPUs, making it a versatile choice for developers looking to implement OpenACC.
- **[Nvidia HPC SDK](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html#gs.zd201n)**: Developed by Nvidia, this compiler is specifically optimized for Nvidia GPUs. It offers robust support for the OpenACC programming model, enabling efficient utilization of GPU resources.
- **[HPE Compiler](https://buy.hpe.com/us/en/software/high-performance-computing-software/high-performance-computing-software/high-performance-computing-software/hpe-cray-programming-environment/p/1012707351)**: Currently, this compiler supports FORTRAN but does not have support for C/C++. It is designed for high-performance computing applications and works well with the OpenACC model.
  

!!! Info "Examples (GNU, Nvidia HPC SDK and HPE): Compilation"

    === "Nvidia HPC SDK"
        ```c
        $ nvc -fast -acc=gpu -gpu=cc80 -Minfo=accel test.c 
        ```

### <u>Questions and Solutions</u>


??? Example "Examples: Hello World"

    === "Serial-version"
        ```c
        //Hello-world-CPU.c	
        #include<stdio.h>
        int main()
        {
          printf("Hello World from CPU!\n");		
          return 0;
        }
        ```


    === "OpenACC-version-parallel"
        ```c
        //Hello-world-parallel.c	
        #include<stdio.h>
        #include<openacc.h>		
        int main()
        { 
        #pragma acc parallel                                                             
          printf("Hello World from GPU!\n");
          return 0;
        }
        ```


    === "OpenACC-version-kernels"
        ```c
        //Hello-world-kernels.c	
        #include<stdio.h>
        #include<openacc.h>		
        int main()
        {
        #pragma acc kernels                            
          printf("Hello World from GPU!\n");
          return 0;
        }
        ```


??? "Compilation and Output"

    === "Serial-version"
        ```c
        // compilation
        $ gcc Hello-world-CPU.c -o Hello-World-CPU
        
        // execution 
        $ ./Hello-World-CPU
        
        // output
        $ Hello World from CPU!
        ```
        
    === "OpenACC-version-parallel"
        ```c
        // compilation
        $ nvc -fast -acc=gpu -gpu=cc80 -Minfo=accel Hello-world-parallel.c -o Hello-World-GPU
        main:
        7, Generating NVIDIA GPU code
        
        // execution
        $ ./Hello-World-GPU
        
        // output
        $ Hello World from GPU!
        ```


    === "OpenACC-version-kernels"
        ```c
        // compilation
        $ nvc -fast -acc=gpu -gpu=cc80 -Minfo=accel Hello-world-kernels.c -o Hello-World-GPU
        main:
        7, Accelerator serial kernel generated
           Generating NVIDIA GPU code
         
        // execution
        $ ./Hello-World-GPU
        
        // output
        $ Hello World from GPU!
        ```


#### <u>Loop</u>

In our second exercise, we will delve into the principles of loop parallelization, a crucial technique in high-performance computing that significantly enhances the efficiency of intensive computations. When dealing with computationally heavy operations within loops, it is often advantageous to parallelize these loops to leverage the full capabilities of multi-core processors or GPUs.

To illustrate this concept, we will begin with a straightforward example: printing **`Hello World from GPU`** multiple times. This will serve as a basis for understanding how to implement loop parallelization effectively.

It is important to note that simply adding directives such as **`#pragma acc parallel`** or **`#pragma acc kernels`** is insufficient for achieving parallel execution of computations. These directives are primarily designed to instruct the compiler to execute the computations on the device, but additional considerations and structure are required to fully exploit parallelization. Understanding how to appropriately organize and optimize loops for parallel execution is essential for maximizing performance in computational tasks. 
    
<figure markdown>
![](../figures/acc-loop.png){align=center width=500}
</figure>
    
    
!!! Info "Loop Constructs"

    === "C/C++"
        ```
        #pragma acc loop [clause-list] new-line
           for loop
        ```
	
    === "FORTRAN"
        ```
        !$acc loop [clause-list]
           do loop
        ```
        

??? "Available clauses for loop"

    === "C/C++ and FORTRAN"
        ```c
        collapse( n )
        gang [( gang-arg-list )]
        worker [( [num:]int-expr )]
        vector [( [length:]int-expr )]
        seq
        independent
        auto
        tile( size-expr-list )
        device_type( device-type-list )
        private( var-list )
        reduction( operator:var-list )
        ```



### <u>Questions and Solutions</u>


??? Example "Examples: Loop (Hello World)"

    === "Serial-version-loop"
        ```c
        //Hello-world-CPU-loop.c	
        #include<stdio.h>
        int main()
        {
          for(int i = 0; i < 5; i++)
            {         
              printf("Hello World from CPU!\n");
            }		
          return 0;
        }
        ```


    === "OpenACC-version-parallel-loop"
        ```c
        //Hello-world-parallel-loop.c	
        #include<stdio.h>
        #include<openacc.h>		
        int main()
        {
        #pragma acc parallel loop
          for(int i = 0; i < 5; i++)
            {                                
              printf("Hello World from GPU!\n");
            }
        return 0;
        }
        ```


    === "OpenACC-version-kernels-loop"
        ```c
        //Hello-world-kernels-loop.c	
        #include<stdio.h>
        #include<openacc.h>		
        int main()
        {
        #pragma acc kernels loop
          for(int i = 0; i < 5; i++)
            {                                
              printf("Hello World from GPU!\n");
            }
        return 0;
        }
        ```


??? "Compilation and Output"

    === "Serial-version-loop"
        ```c
        // compilation
        $ gcc Hello-world-CPU-loop.c -o Hello-World-CPU
        
        // execution 
        $ ./Hello-World-CPU
        
        // output
        $ Hello World from CPU!
        $ Hello World from CPU!
        $ Hello World from CPU!
        $ Hello World from CPU!
        $ Hello World from CPU!                                
        ```
        
        
    === "OpenACC-version-parallel-loop"
        ```c
        // compilation
        $ nvc -fast -acc=gpu -gpu=cc80 -Minfo=accel Hello-world-parallel-loop.c -o Hello-World-GPU
        main:
        5, Generating NVIDIA GPU code
          7, #pragma acc loop gang /* blockIdx.x */
        
        // execution
        $ ./Hello-World-GPU
        
        // output
        $ Hello World from GPU!
        $ Hello World from GPU!
        $ Hello World from GPU!
        $ Hello World from GPU!
        $ Hello World from GPU!                                
        ```


    === "OpenACC-version-kernels-loop"
        ```c
        // compilation
        $ nvc -fast -acc=gpu -gpu=cc80 -Minfo=accel Hello-world-kernels-loop.c -o Hello-World-GPU
        main:
        7, Loop is parallelizable
           Generating NVIDIA GPU code
            7, #pragma acc loop gang, vector(32) /* blockIdx.x threadIdx.x */
        
        // execution
        $ ./Hello-World-GPU
        
        // output
        $ Hello World from GPU!
        $ Hello World from GPU!
        $ Hello World from GPU!
        $ Hello World from GPU!
        $ Hello World from GPU!                                
        ```

