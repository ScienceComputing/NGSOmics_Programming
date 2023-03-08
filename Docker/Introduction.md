# Docker
## What is a Docker?
The Docker is a platform for **consistently** building, running, and shipping applications. The Docker on your development machine can run and function the same way on other test and product machines. The main reasons for using a Docker are: the deployment files are consistent between two machines; the software versions are consistent between two machines; the configuration settings (e.g., ennironment variables) are consistent between two machines. This allows us to consistently run an application on different machines. In a Docker, we can easily package up our application with its dependencies in an isolated environment called `container`. This Docker can be run on any machine. 

## What are the benefits of building containers?
Two applications with different dependencies are packed in isolated environments. They can be installed/deleted or run on the same machine without messing with each other.

## Containers versus virtual machines
A container is an isolated environment for running an application. A virtual machine is an abstraction of a machine (physical hardware). We can run several virtual machines on a real physical machine. For example, we can run the Windows virtual machine and the Linux virtual machine on a Macbook, using the tool `hypervisor`. 

## Hypervisor
The hypervisor is a software to create and manage the virtual machines. For example, VirtualBox, VMware are cross-platform hypervisors which can be run on the Windows, Mac OS, and Linux. The Hyper-v can only be run on the Windows. 

## What are the benefits of building virtual machines?
We can run our applications in isolation.
