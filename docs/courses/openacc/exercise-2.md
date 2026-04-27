#### <u>Data Clauses</u>

Vector addition is a fundamental operation in linear algebra that involves summing two vectors element-wise. Each component of the resulting vector is the sum of the corresponding components from the two input vectors. This example of vector addition highlights two crucial constructs and clauses in OpenACC: compute constructs and data clauses. These include:

- **`#pragma acc parallel loop`**: This directive is useful when the computation involves parallelizing a loop.
- **`#pragma acc kernels loop`**: This directive also applies to loops and enables OpenACC to manage the kernels efficiently.

Data clauses in OpenACC play a pivotal role in seamlessly transferring and managing data between the Central Processing Unit (CPU) and the Graphics Processing Unit (GPU). Below is a description of various data clauses and their usages:

- **`copy`**: Allocates space for a variable on the device, transfers data to the device at the start of the region, copies the data back to the host after the region, and subsequently releases the memory allocated on the device.
- **`copyin`**: Allocates device memory for a variable, transfers data to the device before the region begins, and does not return the data to the host after the region. The device memory is then released. 
- **`copyout`**: Allocates space for a variable on the device but does not copy data to the device before the region. It transfers data back to the host once the region is complete and releases the memory used on the device.
- **`create`**: Allocates memory on the device without transferring any data from the host to the device or vice versa. 
- **`present`**: Indicates that the listed variables already exist on the device, thus requiring no additional action for the transfer. 
- **`deviceptr`**: Utilized for managing data outside of OpenACC, allowing for more controlled access to device pointers.
  
This comprehensive approach to data management enhances the efficiency of GPU computing in high-performance applications.

<figure markdown>
![](../figures/OpenACC-data.png){align=middle, width=750}
<figcaption></figcaption>
</figure>



    
!!! Info "Data Constructs"

    === "C/C++"
        ```
        #pragma acc data [clause-list] new-line
           structured block
        ```
	
    === "FORTRAN"
        ```
        !$acc data [clause-list]
           structured block
        !$acc end data
        ```
        
??? "Available clauses for data"

    === "C/C++ and FORTRAN"
        ```c
        if( condition )
        async [( int-expr )]
        wait [( wait-argument )]
        device_type( device-type-list )
        copy( var-list )
        copyin( [readonly:]var-list )
        copyout( [zero:]var-list )
        create( [zero:]var-list )
        no_create( var-list )
        present(a var-list )
        deviceptr( var-list )
        attach( var-list )
        default( none | present )
        ```



To effectively implement the vector addition example using OpenACC, we need to focus on two specific data clauses. 

1. **Data Transfer for Input Vectors**: The two initialized vectors must be transferred from the host to the device. To achieve this, we will utilize the `copyin` clause, which ensures that the data from the host is available on the device.

2. **Data Transfer for Output Vector**: The product vector, which will store the results of the vector addition, does not require a transfer from the host to the device at the beginning of the computation. However, once the computation is complete, this vector must be transferred back from the device to the host. For this purpose, we will use the `copyout` clause.

By incorporating these data clauses, we can effectively manage the data flow between the host and the device during the execution of the vector addition example. 

In summary, follow these steps to set up the vector addition example with OpenACC: 

1. Use `copyin` to transfer the initialized input vectors to the device.
2. Perform the vector addition on the device.
3. Use `copyout` to transfer the resulting product vector back to the host. 

This approach ensures efficient data handling and optimizes the performance of the application.

Here are the steps for learning the vector addition example:

<figure markdown>
![](../figures/vector_add-external.png){align=middle}
<figcaption></figcaption>
</figure>

 - Allocating the CPU memory for `a`, `b`, and `c` vector
```c
// Initialize the memory on the host
float *restrict a, *restrict b, *restrict c;

// Allocate host memory
a = (float*)malloc(sizeof(float) * N);
b = (float*)malloc(sizeof(float) * N);
c = (float*)malloc(sizeof(float) * N);
```

 - Now, we need to fill in the values for the
    arrays `a` and `b`. 
```c
// Initialize host arrays
for(int i = 0; i < N; i++)
  {
    a[i] = 1.0f;
    b[i] = 2.0f;
  }
```

 - Vector addition kernel function call definition

    ??? "vector addition function call"

        === "Serial-version"
            ```c
            // CPU function that adds two vector 
            void Vector_Addition(float *a, float *b, float *c, int n) 
            {
              for(int i = 0; i < n; i ++)
                {
                  c[i] = a[i] + b[i];
                }
            }
            ```
        
        === "OpenACC-version"
            ```c
            // function that adds two vector 
            void Vector_Addition(float *restrict a, float *restrict b, float *restrict c, int n) 
            {
            #pragma acc kernels loop copyin(a[0:n], b[0:n]) copyout(c[0:n])
            for(int i = 0; i < n; i ++)
              {
                c[i] = a[i] + b[i];
              }
            }	    
            ```




<figure markdown>
![](../figures/vector_add-external-modified.svg) 
<figcaption></figcaption>
</figure>


 - Deallocate the host memory
```c
// Deallocate host memory
free(a); 
free(b); 
free(c);
```

### <u>Questions and Solutions</u>

??? example "Examples: Vector Addition"


    === "Serial-version"
    
        ```c  
        // Vector-addition.c
        
        #include <stdio.h>
        #include <stdlib.h>
        #include <math.h>
        #include <assert.h>
        #include <time.h>
        
        #define N 5120
        #define MAX_ERR 1e-6

        // CPU function that adds two vector 
        float * Vector_Addition(float *a, float *b, float *c, int n) 
        {
          for(int i = 0; i < n; i ++)
            {
              c[i] = a[i] + b[i];
            }
          return c;
        }

        int main()
        {
          // Initialize the memory on the host
          float *a, *b, *c;       
  
          // Allocate host memory
          a = (float*)malloc(sizeof(float) * N);
          b = (float*)malloc(sizeof(float) * N);
          c = (float*)malloc(sizeof(float) * N);
  
          // Initialize host arrays
          for(int i = 0; i < N; i++)
            {
              a[i] = 1.0f;
              b[i] = 2.0f;
            }
    
          // Start measuring time
          clock_t start = clock();

          // Executing CPU function 
          Vector_Addition(a, b, c, N);

          // Stop measuring time and calculate the elapsed time
          clock_t end = clock();
          double elapsed = (double)(end - start)/CLOCKS_PER_SEC;
        
          printf("Time measured: %.3f seconds.\n", elapsed);
  
          // Verification
          for(int i = 0; i < N; i++)
            {
              assert(fabs(c[i] - a[i] - b[i]) < MAX_ERR);
            }
          printf("PASSED\n");
    
          // Deallocate host memory
          free(a); 
          free(b); 
          free(c);
   
          return 0;
        }
        ```

    === "OpenACC-template"
    
        ```c  
        // Vector-addition-template.c
	
        #include <stdio.h>
        #include <stdlib.h>
        #include <math.h>
        #include <assert.h>
        #include <time.h>
        #include <openacc.h>	


        #define N 5120
        #define MAX_ERR 1e-6


        // GPU function that adds two vectors 
        // function that adds two vector 
        void Vector_Addition(float *restrict a, float *restrict b, float *restrict c, int n) 
        {

        // add here either parallel or kernel plus data map clauses
        #pragma acc 
        for(int i = 0; i < n; i ++)
           {
             c[i] = a[i] + b[i];
           }
        }

        int main()
        {
          // Initialize the memory on the host
          float *restrict a, *restrict b, *restrict c;       
  
          // Allocate host memory
          a = (float*)malloc(sizeof(float) * N);
          b = (float*)malloc(sizeof(float) * N);
          c = (float*)malloc(sizeof(float) * N);
  
          // Initialize host arrays
          for(int i = 0; i < N; i++)
            {
              a[i] = 1.0f;
              b[i] = 2.0f;
            }
    
          // Start measuring time
          clock_t start = clock();

          // Executing CPU function 
          Vector_Addition(a, b, c, N);

          // Stop measuring time and calculate the elapsed time
          clock_t end = clock();
          double elapsed = (double)(end - start)/CLOCKS_PER_SEC;
        
          printf("Time measured: %.3f seconds.\n", elapsed);
  
          // Verification
          for(int i = 0; i < N; i++)
            {
              assert(fabs(c[i] - a[i] - b[i]) < MAX_ERR);
            }

          printf("PASSED\n");
    
          // Deallocate host memory
          free(a); 
          free(b); 
          free(c);
   
          return 0;
        }
        ```
	
    === "OpenACC-version"
    
        ```c  
        // Vector-addition-openacc.c
        
        #include <stdio.h>
        #include <stdlib.h>
        #include <math.h>
        #include <assert.h>
        #include <time.h>
        #include <openacc.h>

        #define N 5120
        #define MAX_ERR 1e-6


        // function that adds two vector 
        void Vector_Addition(float *restrict a, float *restrict b, float *restrict c, int n) 
        {
        // or #pragma acc kernels loop copyin(a[0:n], b[0:n]) copyout(c[0:n])
        #pragma acc kernels loop copyin(a[0:n], b[0:n]) copyout(c[0:n])
        for(int i = 0; i < n; i ++)
           {
            c[i] = a[i] + b[i];
           }
        }

        int main()
        {
          // Initialize the memory on the host
          float *restrict a, *restrict b, *restrict c;       
  
          // Allocate host memory
          a = (float*)malloc(sizeof(float) * N);
          b = (float*)malloc(sizeof(float) * N);
          c = (float*)malloc(sizeof(float) * N);
  
          // Initialize host arrays
          for(int i = 0; i < N; i++)
            {
              a[i] = 1.0f;
              b[i] = 2.0f;
            }
    
          // Start measuring time
          clock_t start = clock();

          // Executing CPU function 
          Vector_Addition(a, b, c, N);

          // Stop measuring time and calculate the elapsed time
          clock_t end = clock();
          double elapsed = (double)(end - start)/CLOCKS_PER_SEC;
          
          printf("Time measured: %.3f seconds.\n", elapsed);
  
          // Verification
          for(int i = 0; i < N; i++)
            {
              assert(fabs(c[i] - a[i] - b[i]) < MAX_ERR);
            }

          printf("PASSED\n");
    
          // Deallocate host memory
          free(a); 
          free(b); 
          free(c);
   
          return 0;
        }
        ```


??? "Compilation and Output"

    === "Serial-version"
        ```c
        // compilation
        $ gcc Vector-addition.c -o Vector-Addition-CPU
        
        // execution 
        $ ./Vector-Addition-CPU
        
        // output
        $ ./Vector-addition-CPU 
        PASSED
        ```
        
    === "OpenACC-version"
        ```c
        // compilation
        $ nvc -fast -acc=gpu -gpu=cc80 -Minfo=accel Vector-addition-openacc.c -o Vector-Addition-GPU
        Vector_Addition:
        12, Generating copyin(a[:n]) [if not already present]
            Generating copyout(c[:n]) [if not already present]
            Generating copyin(b[:n]) [if not already present]
        14, Loop is parallelizable
            Generating NVIDIA GPU code
            14, #pragma acc loop gang, vector(128) /* blockIdx.x threadIdx.x */

        // execution
        $ ./Vector-Addition-GPU
        
        // output
        $ ./Vector-addition-GPU
        PASSED
        ```


??? Question "Question"

    - Please try other data clauses for	different applications and get familiarised with them.

