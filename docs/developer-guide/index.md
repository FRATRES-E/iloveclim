# Developer guide

## Before you start

Contribution to the iLOVECLIM codebase requires special rights. FRATRES is built with three levels of permissions:

1. The code itself is public access with an Apache 2.0 licensing. Anyone can evaluate the code in a spirit of scientific openess. 

2. Accessing the full build system requires access to the iloveclim-users team, essentially giving full read access to the i-operator submodule. 

3. Contributing code to iLOVECLIM can be done either through a clone, pull request and review by the teams of developpers or, if member of the iloveclim-dev team, by a direct push to the main branch. Either way, contact is strongly advised with the members of the iloveclim-dev to ensure useability and non-duplication. 

## Outline

The following is structured as follows:

- The general architecture of the model is described in [code architecture](architecture.md)

- Coding convention and coding templates are then [explained and provided](conventions.md)

- Being a github based project, the [github workflow](git-workflow.md) is the place to go when you are starting to develop a new feature. A few conventions for commits are provided there in addition as well.

- At a greater abstraction level, [component addition](adding-a-component.md) might become a topic of interest. This comes with higher level of complexity into the coupling parts and the software build workflow. 

- Finally some recommendation (incomplete) are provided so as to provide [testing and CI](testing.md) links. These will probably grow over time. 

 
