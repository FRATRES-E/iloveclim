# Obtaining the code of iLOVECLIM

## Introduction

The code of iLOVECLIM consists of two separate repository: the main public iLOVECLIM repository in FRATRES-E, a public repository, where anyone can always check the code, propose changes through pull requests etc. The second one is the i-operator which is the machinery to get the model up and running. This section of the code is restricted to the active members of the iLOVECLIM group and is not public, you need to be a member of the proper group to check the code out. 

The *i-operator* is an independent github repository which can be checked out as a submodule in iLOVECLIM (or sub git repository if you wish). In the following, the installation guide assumes that you have the right permissions to do so. 

!!! warning "Github tokens"
    For checking out the private *i-operator*, github will ask you what looks like a username and a password. Though not specified explicitly on the command line, what is expected is a **token** that you are expected to have generated in your account and kept private.



## iLOVECLIM checkout

Clone the git iloveclim directory to a fresh installation (replace "your_new_iloveclim" by your own directory name)

    git clone --recursive https://github.com/FRATRES-E/iloveclim.git your_new_iloveclim 

Just a comment : the iloveclim code is now public in FRATRES-E, but all our beloved installation and running tools are not. Hence you need to be part of the iloveclim-users team to get the submodule "i-operator" that does that part of the work.

Get the standard input data (not part of iLOVECLIM anymore) at LSCE, something like, within the iLOVECLIM directory:

    wget https://dods.lsce.ipsl.fr/iloveclim-dmr/ilcm-inps/inputs_iloveclim-2025-10-08.tar.bz2

Unpack that previous package in the iLOVECLIM main directory with:

    tar -xvjf inputs_iloveclim-2025-10-08.tar.bz2

You are now ready for the install process.

## iLOVECLIM Installation

First, we need to generate the makefile and make.macros/libinc files through lipas. Go to your LIPaS directory and perform something like:

    ./lipas.sh -g absolute_path_towards_iloveclim_dir/i-operator/i-operator.toml

Then get back to the iLOVECLIM directory, in the i-operator subfolder. 
Create an auto.lipas file using:

    ls ${HOME}/.lipas/config/*/*/conf.* > config/auto.lipas

(SANITY CHECK: auto.lipas should contain a path to one existing file)

Then, back to your iLOVECLIM folder, perform:

    ./i-operator/Niloveclim_install-FRATRES.sh

... that's it!!!
