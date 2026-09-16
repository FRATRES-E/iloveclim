# Base model update

As clearly apparent from the installation page itself, the code of the model is under a github repository and is under active development. This is not a polished finalized software but a research tool. It is thus probable that you encounter bugs that need to be fixed and that you will need to update your installation. We have seen that the model is split in two: a public repository and a private submodule. Strategies to update both are necessarily different. 

## Update to the iLOVECLIM code itself

If the update you seek to include in your iLOVECLIM installation is only affecting the source code of the model (and not the install itself) then it is rather simple: go to your iLOVECLIM installation directory and update the code with a:

```bash
git pull
```

that's it. Or that would it if you would not have some local modifications of the code that created some conflicts, but that, as commonly accepted, is another story (see page to be written on basic git usage).

## Update to the i-operator

Often enough, you will need to update the model and its installer / run script. For this you need to coherently update both the iLOVECLIM and its `i-operator`. As the latter is a git submodule, you need the slightly more complex command:

```bash
git pull && git pull --recurse-submodules && git submodule update
```

You will be prompted for a user name and password (which in fact is your beloved _git token_) and you will be done with the checkout. You will still need to perform an install update of the model, that is re-running the:
```bash
./i-operator/Niloveclim_install-FRATRES.sh
```
This will generate a new version of the `run-me.sh` and you will be all set.
