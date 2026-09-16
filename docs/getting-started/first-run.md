# Your first run

!!! warning "Version"
    At time of writing (2026-09-16) the run system is what will be in the future the *legacy* system. Ongoing work is being done to render all this section more flexible. **Stay tuned!!**

## First things first: understanding the run command

Since iLOVECLIM is based on a number of components and we like to keep things along the [KISS principle](https://en.wikipedia.org/wiki/KISS_principle), there is [One Script to Rule them All](https://tolkiengateway.net/wiki/Ring-inscription) that is full of options and is doing (most or all) the work for you. 

### Basic use of the run script
If you successfully finalized the installation, you should have a `run-me.sh` in your iLOVECLIM base directory. This is a symbolic link to the file `i-operator/bin/run-iloveclim.gen` that has been generated during the installation process. 

To test first if things are allright, the simplest possible use is:
```bash
./run-me.sh
```
This will perform a standard run, with standard options as is set in the source code by default; it will do a coupled run with atmosphere--ocean--vegetation activated (land-surface is implicitly also there) at Pre-Industrial conditions for a length of one year, starting from a previous equilibrium. This is the standard run that should always work, if not you should shout it out loud to the iLOVECLIM mailing list! (see page **Troubleshooting**). 

### Analyzing the options
If you are curious about what the script did for you in the background (as you should :smile:) you can run the following to list the options:
```bash
./run-me.sh -h
```
to get the online help. You should obtain something like:
```

	 	 	 	 iLOVECLIM 
	 	 	 	 (http://forge.ipsl.jussieu.fr/ludus) 	 

[MAIN] iLOVECLIM CALL: ./run-me.sh -h
Usage: ./run-me.sh OPTIONS
OPTIONS:
    -h                     Print this help message (+Queue info)
    -l <run_label>         A label identifying the experiment                     (Default: test)
    -s <start_year>        The year in which the run starts (B.P.> celest)        (Default: 3000)
    -n <num_years>         The length of the run in model years                   (Default: 1)
    -r <restart_interval>  The number of years between writing a restart file     (Default: 1)
    -f <scenario_label>    A label identifying the scenario                       (Default: default)
    -z                     Setup run, but do not run                              (Default: 0)
    -e                     Model sensitivity parameter version                    (Default: E00)
    -v <level>             verbose 1: basic info, 2: allinfo                      (Default: 1)

NOT DEFAULT OPTIONS:
    -k                     remove wkdir & data dir when they exist                (Default: no)
    -t <scratch_dir>       run on a local/scratch dir, retrieve results at end    (Default: none)

 ===================== Additional options for non-standard runs ================================
 ====== Files provided in this part are saved in /data/<scenario_label>/sub-directories  =======
 ===============================================================================================

    -P <directory path>   A directory containing modified parameter files        (Default: None)
    -F <directory path>   A directory containing modified source files           (Default: None)
    -I <directory path>   A directory containing modified *input*data files      (Default: None)
    -S <directory path>   A directory containing non-standard *re-Start* files   (Default: None)
```
The website width does not allow to display it, but there are Default values for each value in this list. The standard default is equivalent to `-l test -n 1 -r 1` etc. 

Would you wish to run a longer experiment in the same test context you could for example run:
```bash
./run-me.sh -l test_long -n 100 -r 100 
```
to run a hundred years of equilibrium on your machine. 

!!! tip "Restart intervals"
    The model will write data to allow performing a restart at certain indicated points. By default, the value that controls it (`-r 1`) is set to one. In the last example above, I have set this value to `100`, which means that the writing of the restart will be done at the end of the 100 years run. You can set whatever you wish, but if the restart interval is greater than the run length you will have no restart written and if you do it too frequently, it will take quite some space on disk. 

### Additional options
In the previous options lookout, you can see that there are four different options that allow replacing on the fly some standard files in a model run:

=== "`-P` Parameter files"
    This allows using a non standard parameter file. The argument is a directory that should contain parameter files with standard names but customized content. They will be copied in the _run_ directory and used in place of the default one with the same name (simple overwriting). 

    
=== "`-F` _FORTRAN_ (source) files"
    This allows using a non standard code source files (_FORTRAN_ language, hence `-F`). The argument is a directory that should contain source files with standard names but customized content. They will be copied in the _compil_ directory and used in place of the default one with the same name (simple overwriting). They will be copied before the compilation process. In the case where you put a source file that is non standard, it will be also copied (simple globbing of _FORTRAN_ extensions and since the model does not have a unique list of source files, it will be including in the compilation pattern through automatic dependency solving.

=== "`-I` Inpudtdata files"
    This allows using a non standard inputdata file. The argument is a directory that should contain inputdata files with standard names but customized content. They will be copied in the _run/inputdata_ directory and used in place of the default one with the same name (simple overwriting). **Warning:** the inputs have a standard directory structure by model sub-components, you need to have a mirror structure in that additional directory in order to ensure proper operation. 

=== "`-R` Restartdata state"
    This allows using a non standard restartdata directory. The argument is a directory that should contain a full restart state. This is maybe the option that is used most. The restart state will be copied in the _run_ directory and used as is. **Warning:** this will only work if you also provide the `-s` option to the script with a restart number that matches your restart state (for restart state `res001234`, you need to add `-s 1234`, without leading zeros. 


## Looking at your first results

### In the event of a nice run

Once your run is happily completed, the model has generated a number of output files that you might want to have a look at. They are located either in the `data/name_of_your_run/output00XXXX` or in `wkdir/name_of_your_run/run/outputdata` depending on your exact setup. On the long run, we are moving towards generalization of the second solution. 

The `output` directory (whatever its name) contains a series of subdirectories with the output of the different components of the model. The simple default run yields:

```
atmos   climate_indices  downscaling  icebergs  land   permafrost  vegetation
carbon  coupler          globals      ism       ocean  sediments
```

Each of these subdirectories contains the output of the component named (provided it was activated, if not it is empty) and the `globals` provide some generic diagnostics of the whole model. Most of the output is in netCDF format and on the internal generic grid of the component (for the ocean for example it is on a pole rotated grid that is a little harder to plot -- but no panic!). 

### In the event of a nasty run

If you run directly on the front machine (which, should I stress it, you should **not** do), your output will still be in the given `run` directory with incomplete (or not at all output). If you are running on a cluster with a job scheduler, please refer to the dedicated page for more information ; in principle the setup of the model is such that it tries to retrieve what existed at the time of crash - if at all.
