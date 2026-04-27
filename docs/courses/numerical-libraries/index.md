# Introduction to numerical libraries

The Luxembourg SuperComputing Competence Center will host a half-day online course on numerical libraries. The course will be divided into two parts: the first will cover theoretical concepts, while the second will focus on practical, hands-on challenges using the MeluXina supercomputer.

### Who should attend
This course is ideal for current and prospective users of large hybrid CPU/GPU clusters and supercomputers who develop application that rely on numerical algorithms such as the vector and matrix operations implemented in BLAS and the basic linear algebra operations implemented in LAPACK.

### What will you learn and how
Participants in this course will learn how to use numerical libraries in C/C++ programs. Calling library function and linking with library object files will be covered first, followed by a detailed investigation of how the effects of caching affect the performance of libraries highly optimized to exploit the caching.

A presentation will cover the data structures and the algorithms used in BLAS and LAPACK, and how they are designed to exploit the cache. Accelerated implementations such as cuBLAS and MAGMA will also be presented, including an analysis of data movement to and from the accelerator that tends to determine the performance of accelerated libraries. Then, the participants will have the opportunity to play with demonstrative examples using BLAS and LAPACK in the practical session. Examples with accelerated libraries will also be available utilizing the GPU partition of MeluXina.

### Learning Outcomes
####By the end of the course, participants will be able to:
 - Understand the design of numerical libraries, including:
    - Matrix data structures for dense and sparse matrices (CSR, ELL, COO)
    - Algorithmic implementation of linear algebra operations and how it takes advantage of the cache
 - Compile efficient implementations of numerical libraries, including:
    - Enabling implementation specific optimizations
    - Inspecting object files to determine which components of the library interface they support
    - Using libraries in compilations with build automation systems such as CMake
 - Efficiently exploit cache to leverage the best performance, with methods such as:
    - Cache aware programming
    - Cache alignment
 - Efficient use libraries optimized for accelerators, including:
    - Understanding the difference between cache and accelerator-optimized libraries
    - Efficiently manage data movement to and from the accelerator

### Prerequisites

Priority will be given to participants with solid experience in C/C++ and/or FORTRAN. Some knowledge about the steps involved in compilation and linking would be helpful but not necessary.

### Computing Resources

Participants will have access to the MeluXina supercomputer CPU and GPU partitions during the training. For more information about MeluXina, please refer to the MeluXina Overview and the MeluXina – Getting Started Guide.

### Agenda

This half-day course will be conducted online in Central European Time (CET) on February 14th, 2025 (9:00 AM – 13:00 PM CET).

### Schedule: 
 - 9:00 AM – 9:50 AM: Lecture Part 1: The data structure and algorithms of BLAS and LAPACK
 - 9:50 AM – 10:00 AM: Break
 - 10:00 AM – 10:50 AM: Lecture Part 2: Practical aspects of working with libraries
 - 10:50 AM – 11:00 AM: Break
 - 11:00 AM – 11:50 AM: Lecture Part 3: Effects of cache and accelerated implementations
 - 11:50 AM – 12:00 PM: Break
 - 12:00 PM – 12:50 PM: Practical session
 - 12:50 PM – 13:00 PM: Q&A

## Important: Limited spots are available!

Contact people for more info: <br />
Georgios KAFANAS, georgios.kafanas@uni.lu 
