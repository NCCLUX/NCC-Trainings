Profiling is an important task to be considered when a computer code is written. Writing parallel code is less challenging, but making it more efficient on a given parallel architecture is challenging. Moreover,  from the programming and programmer’s perspective, we want to know where the code spends most of its time. In particular, we would like to know if the code (given algorithm) is compute bound, memory bound, cache misses, memory leak, proper vectorisation, cache misses, register spilling, or hot spot (time-consuming part in the code). Plenty of tools are available to profile a scientific code (computer code for doing arithmetic computing using processors). However, we will focus on a few of the widely used tools.

 - [AMD uProf](https://www.amd.com/content/dam/amd/en/documents/developer/uprof-v4.0-gaGA-user-guide.pdf)
 - [ARM Forge](https://developer.arm.com/documentation/101136/22-1-3/Performance-Reports?lang=en)
 - [Intel tools](https://www.intel.com/content/www/us/en/developer/tools/oneapi/vtune-profiler.html)

####<u>[ARM Forge](https://developer.arm.com/documentation/101136/22-1-3/Performance-Reports?lang=en)</u>
Arm Forge is another standard commercial tool for debugging, profiling, and analysing scientific code on the massively parallel computer architecture. They have a separate toolset for each category with the common environment: DDT for debugging, MAP for profiling, and performance reports for analysis. It also supports the MPI, UPC, CUDA, and OpenMP programming models for different architectures with a variety of compilers. DDT and MAP will launch the GUI, where we can interactively debug and profile the code. Meanwhile, `perf-report` will provide the analysis results in `.html` and `.txt` files.


??? Info "Example: ARM Forge"

    === "C/C++"
        ```c
        # compilation with debugging tool
        $ gcc test.c -g -fopenmp
        # execute and profile the code
        $ map --profile --no-mpi ./a.out
        # open the profiled result in GUI
        $ map xyz.map
        
        # for debugging
        $ ddt ./a .out
        
        # for profiling
        $ map ./a .out
        
        # for analysis
        $ perf-report ./a .out
        ```

    === "FORTRAN"
    	```c
        # compilation 
        $ gfortran test.f90 -fopenmp
        # execute and profile the code
        $ map --profile --no-mpi ./a.out
        # open the profiled result in GUI
        $ map xyz.map
        
        # for debugging
        $ ddt ./a .out
        
        # for profiling
        $ map ./a .out

        # for analysis
        $ perf-report ./a .out
        ```
    <figure markdown>
    ![](../figures/Arm-Forge.png){align=center}
    <figcaption></figcaption>
    </figure>

####<u>[Intel tools](https://www.intel.com/content/www/us/en/developer/tools/oneapi/vtune-profiler.html)</u>

##### Intel Application Snapshot
Intel Application Performance Snapshot tool helps to find essential performance factors and the metrics of CPU utilisation, memory access efficiency, and vectorisation.
`aps -help` will list out profiling metrics options in APS
     
<figure markdown>
![](../figures/APS_OpenMP_flow_chart.png){align=center}
<figcaption></figcaption>
</figure>

??? Info "Example: APS"

    === "C/C++"
        ```c
        # compilation
        $ icc -qopenmp test.c
        
        # code execution
        $ aps --collection-mode=all -r report_output ./a.out
        $ aps-report -g report_output                        # create a .html file
        $ firefox report_output_<postfix>.html               # APS GUI in a browser
        $ aps-report report_output                           # command line output
        ```

    === "FORTRAN"
    	```c
        # compilation
        $ ifort -qopenmp test.f90
        
        # code execution
        $ aps --collection-mode=all -r report_output ./a.out
        $ aps-report -g report_output                        # create a .html file
        $ firefox report_output_<postfix>.html               # APS GUI in a browser
        $ aps-report report_output                           # command line output
        ```

    <figure markdown>
    ![](../figures/OpenMP_APS.png){align=center}
    <figcaption></figcaption>
    </figure>

##### Intel Inspector

Intel Inspector detects and locates the memory, deadlocks, and data races in the code.
For example, memory access and memory leaks can be found.

??? Info "Example: Intel Inspector"

    === "C/C++"
        ```c
        # compile the code
        $ icc -qopenmp example.c
        # execute and profile the code
        $ inspxe-cl -collect mi1 -result-dir mi1 -- ./a.out
        $ cat inspxe-cl.txt
        # open the file to see if there is any memory leak
        === Start: [2020/12/12 01:19:59] ===
        0 new problem(s) found
        === End: [2020/12/12 01:20:25] ===
        ```

    === "FORTRAN"
    	```c
        # compile the code
        $ ifort -qopenmp test.f90
        # execute and profile the code
        $ inspxe-cl -collect mi1 -result-dir mi1 -- ./a.out
        $ cat inspxe-cl.txt
        # open the file to see if there is any memory leak
        === Start: [2023/05/10 01:19:59] ===
        0 new problem(s) found
        === End: [2020/05/10 01:20:25] ===
        ```


##### Intel Advisor

Intel Advisor: a set of collection tools for the metrics and traces that can be used for further
tuning in the code. `survey`: analyse and explore an idea about where to add efficient vectorisation.


??? Info "Example: Intel Advisor"

    === "C/C++"
        ```c
        # compile the code
        $ icc -qopenmp test.c
        # collect the survey metrics
        $ advixe-cl -collect survey -project-dir result -- ./a.out
        # collect the report
        $ advixe-cl -report survey -project-dir result
        # open the gui for report visualization
        $ advixe-gui
        ```

    === "FORTRAN"
        ```c
        # compile the code
        $ ifort -qopenmp test.90
        # collect the survey metrics
        $ advixe-cl -collect survey -project-dir result -- ./a.out
        # collect the report
        $ advixe-cl -report survey -project-dir result
        # open the gui for report visualization
        $ advixe-gui
        ```

    <figure markdown>
    ![](../figures/Advisor.png){align=center}
    <figcaption></figcaption>
    </figure>
    

##### Intel VTune

 - Identify the time-consuming part of the code.
 - Also, identify cache misses and latency.


??? Info "Example: Intel VTune"

    === "C/C++"
        ```c
        # compile the code
        $ icc -qopenmp test.c
        # execute the code and collect the hotspots
        $ amplxe-cl -collect hotspots -r amplifier_result ./a.out
        $ amplxe-gui
        # open the GUI of the VTune amplifier
        ```

    === "FORTRAN"
        ```c
        # compile the code
        $ ifort -qopenmp test.90
        # execute the code and collect the hotspots
        $ amplxe-cl -collect hotspots -r amplifier_result ./a.out
        $ amplxe-gui
        # open the GUI of the VTune amplifier
        ```
	
    <figure markdown>
    ![](../figures/Vtune.png){align=center}
    <figcaption></figcaption>
    </figure>
    
    `amplxe-cl` will list out the analysis types and `amplxe-cl -hlep` report will list out available reports in VTune.


####<u>[AMD uProf](https://www.amd.com/content/dam/amd/en/documents/developer/uprof-v4.0-gaGA-user-guide.pdf)</u>
AMD uProf profiler follows a statistical sampling-based approach to collect profile data to identify
the performance bottlenecks in the application.

??? Info "Example: AMD uProf"

    === "C/C++"
        ```c
        # compile the code
        $ clang -fopenmp test.c
        $ AMDuProfCLI collect --trace openmp --config tbp --output-dir solution ./a.out -d 1
        ```

    === "FORTRAN"
        ```c
        # compile the code
        $ flang -fopenmp test.90
        $ AMDuProfCLI collect --trace openmp --config tbp --output-dir solution ./a.out -d 1
        ```
