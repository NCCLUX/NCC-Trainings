# Environment management and containers for HPC systems

Environment management and containerization in shared computing resources like HPC systems and workstations allow the installation and execution of programs in isolation of other system and user applications. Environment management and containerization encompass a collection of techniques base on features of the operating system and the kernel which are used for a diverse set of tasks such as

- sourcing and installing software
- managing conflicting dependencies of software, and
- creating reproducible environments.

This course will start by providing and overview of the operating systems and kernel features that are used in environments and containers. Then, we cover the basic skills for creating environments and containers so that you can develop and deploy your software in shared systems. Finally, we cover the necessary techniques that allow the creation of reproducible environments and reproducible binaries for containers. Reproducibility is particularly important in scientific settings, where workflows for simulations and results analysis can be encapsulated and distributed in a single container or container recipe to ensure the reproducibility and traceability of results.

## Target audience

The course targets users of HPC systems that need to develop or deploy applications in a consistent manner across multiple systems.

## Agenda

- Overview of the kernel and operating system features used to manage software environments.
- Basic software environments for software installation and dependency management.
- Introduction to container for the creation and distribution of consistent software environments.
- Advanced management of system resources in containers: accessing storage, networks, and accelerators.
- Creating reproducible environments and containers.

## Requirements

- Having an HPC account to access the cluster.
- Basic shell usage skills.
- Basic scripting skills and the ability to use a text editor in the cluster.
