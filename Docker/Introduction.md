# Docker
## What is a Docker?
The Docker is a platform for **consistently** building, running, and shipping applications. The Docker on your development machine can run and function the same way on other test and product machines. The main reasons for using a Docker are: the deployment files are consistent between two machines; the software versions are consistent between two machines; the configuration settings (e.g., ennironment variables) are consistent between two machines. This allows us to consistently run an application on different machines. In a Docker, we can easily package up our application with its dependencies in an isolated environment called `container`. This Docker can be run on any machine. 

## What are the benefits of building containers?
Two applications with different dependencies are packed in isolated environments. They can be installed/deleted or run on the same machine without messing with each other.

## Containers versus virtual machines
A container is an isolated environment for running an application. It is more lightweight and does not require the full operating system, compared to a virtual machine. All containers on a single machine shares the operating system of the host. A container can start more quickly than a virtual machine, because the operating system has already started on the host. A container also needs fewer hardware resources (e.g., number of CPU cores, memory, disk space) on the host than a virtual machine. On a single host, we can run tens or even hundres of containers, side by side. A virtual machine is an abstraction of a machine (physical hardware). We can run several virtual machines on a real physical machine. For example, we can run the Windows virtual machine and the Linux virtual machine on a Macbook, using the tool `hypervisor`. 

## Hypervisor
The hypervisor is a software to create and manage the virtual machines. For example, VirtualBox, VMware are cross-platform hypervisors which can be run on the Windows, Mac OS, and Linux. The Hyper-v can only be run on the Windows. 

## What are the benefits of building virtual machines?
We can run our applications in isolation. On the same physical machine, we can have two different virtual machines, each of them running a different application with different dependencies it needs. 

## What are the problems of building virtual machines?
Each virtual machine needs a full copy of operating system that needs to be licensed, patched, and monitored. Therefore, virtual machines are slow to start, since the entire operating system needs to be loaded. Each virtual machine is resource intensive because it takes the slice of the actual hardware resources such as CPU, memory, and disk space. For example, if the physical machine has a memory of 16 gigabytes, this memory has to be divided between different virtual machines. We can determine how much memory to allocate to each virtual machine. However, as the virtual machine is resource intensive, we are limited by the number of virtual machines to run on a physical machine.
