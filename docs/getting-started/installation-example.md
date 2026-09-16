## iLOVECLIM Installation

First, we need to generate the makefile and make.macros/libinc files through lipas. Go to your LIPaS directory and perform something like:
```bash
./lipas.sh -g absolute_path_towards_iloveclim_dir/i-operator/i-operator.toml
```
Then get back to the iLOVECLIM directory, in the i-operator subfolder. 
Create an auto.lipas file using:
```bash
ls ${HOME}/.lipas/config/*/*/conf.* > config/auto.lipas
```

(SANITY CHECK: auto.lipas should contain a path to one existing file)

Then, back to your iLOVECLIM folder, perform:

```bash
./i-operator/Niloveclim_install-FRATRES.sh
```

... that's it!!!

