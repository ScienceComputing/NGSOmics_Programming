# Docker
## What is a Docker?
The Docker is a platform for **consistently** building, running, and shipping applications. The Docker on your development machine can run and function the same way on other test and product machines. The main reasons for using a Docker are: the deployment files are consistent between two machines; the software versions are consistent between two machines; the configuration settings (e.g., ennironment variables) are consistent between two machines. This allows us to consistently run an application on different machines. In a Docker, we can easily package up our application with its dependencies in an isolated environment called `container`. This Docker can be run on any machine. 

## What are the benefits of building containers?
Two applications with different dependencies are packed in isolated environments. They can be installed/deleted or run on the same machine without messing with each other. If we are done with one application, we can safely remove it with all its dependencies in one go, cleaning up our machine. We do not need to worry if the removal of these dependencies will impact the running of other applications. 

## Containers versus virtual machines
A container is an isolated environment for running an application. It is more lightweight and does not require the full operating system, compared to a virtual machine. All containers on a single machine `shares the operating system` of the host, or more accurately, all containers on a single machine `shares the kernel` of the host. This means we only need to license, patch, and monitor a single operating system. A container can start more quickly than a virtual machine, because the operating system has already started on the host. A container also needs fewer hardware resources (e.g., number of CPU cores, memory, disk space) on the host than a virtual machine. On a single host, we can run tens or even hundres of containers, side by side. A virtual machine, VM, is an abstraction of a machine (physical hardware), which contains a full-blown operating system. We can run several virtual machines on a real physical machine. For example, we can run the Windows virtual machine and the Linux virtual machine on a Macbook, using the tool `hypervisor`. 

## Kernel
Kernel is the core of an operating system. It manages all applications and hardware resources such as memory and CPU. Every operating system has its own kernel. These kernels have different APIs. That's why we cannot run a Windows application on a Linux operating system, because under the hood, this application cannot talk to the Linux kernel. On a Linux machine, we can only run Linux containers, as Linux containers share the Linux kernel. At present, a Windows machine is built with both the Windows kernel and Linux kernel, so on this machine, we can run both Windows containers and Linux containers. Our windows containers share the windows kernel and our linux containers share the linux kernel. On a Linux machine, we can only run Linux containers. On a Mac machine, we build a lightweight Linux virtual machine to run Linux containers such as Docker. 

## Hypervisor
The hypervisor is a software to create and manage the virtual machines. For example, VirtualBox and VMware are cross-platform hypervisors which can be run on the Windows, Mac OS, and Linux. The Hyper-v can only be run on the Windows. 

## What are the benefits of building virtual machines?
We can run our applications in isolation inside the virtual machine. On the same physical machine, we can have two different virtual machines, each of them running a different application with different dependencies it needs. 

## What are the problems of building virtual machines?
Each virtual machine needs a full copy of operating system that needs to be licensed, patched, and monitored. Therefore, virtual machines are slow to start, since the entire operating system needs to be loaded. Each virtual machine is resource intensive because it takes the slice of the actual physical  hardware resources such as CPU, memory, and disk space. For example, if the physical machine has a memory of 16 gigabytes, this memory has to be divided between different virtual machines. We can determine how much memory to allocate to each virtual machine. However, as the virtual machine is resource intensive, we are limited by the number of virtual machines to run on a physical machine.

## Docker architecture
The Docker uses the client-server architecture where it has a client component that converses with a server component using a rest API. The server is also called "Docker Engine". It sits on the background and takes care of building and running Docker containers. Techinically, **a container is a process, but it's a special kind of process because it has its own file system which is provided by the image**. So our application gets loaded inside a container or a process and this is how we run our application locally on our development machine. Rather than launching the application directly and executing it within a conventional process, we instruct Docker to execute it within a container, thereby creating an isolated environment.

## Install Docker
https://docs.docker.com/get-docker/

On mac and windows, we have docker desktop which is the combination of docker engine plus other tools. 

## Check Docker version
`docker version`

## Docker file
A docker file is a plain text file which includes instructions that docker employs to bundle this application into an image. This image contains everything an application rests upon to run, for example, an operating system, a runtime environment like node or python, application files, third-party libraries and environment variables. 

With docker we no longer need to maintain long complex release documents that have to be precisely followed, all the instructions for building an image of an application are written in a docker file. With that we can package up our application into an image and run it virtually anywhere. 

## Run a Docker container
`docker run ...`

## Docker Hub
This is a registry where docker images can be pushed into. The Docker Hub to docker is like the GitHub to git. It's a storage for docker images that anyone can use so once our application image is on Docker Hub then we can download it on any machines. 
