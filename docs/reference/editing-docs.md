# A group effort in being: editing the docs

The docs of the iLOVECLIM group are part of the iLOVECLIM git repository. These are simply a collection of markdown files that are hosted in the `docs` repository of the model. Anyone with developper access can edit the docs and commit changes. Anybody else can create a pull request to modify the docs accordingly. 

!!!warning Developper mode only
    In the following, I assume that you have the developper rights on the iLOVECLIM `git` repository. Contact the repository owners if you want to know more. 

## Editing the docs locally

When you have a local copy of the iLOVECLIM model, you have a copy of the [docs.iloveclim.eu](http://docs.iloveclim.eu) website with you. The `docs` repository contains everything in the form of a structure with subdirectories containing markdown files. If you want to update some webpage, simply edit the text content with your favorite editor (if you want a markdown editor, have a look at [marktext](https://github.com/marktext/marktext) for example). Once you are happy with your content, you can add the files to a `git commit` phase and then make a `git push` to update the github website. From there on, _github Actions_ will take over to generate a new version of the website and make it available in the proper place. 

!!!warning "A note of caution" 
    Whatever you push in the github repository will be directly deployed on the dosc website publicly. Act accordingly but wisely.

## Testing the docs on a local version

Most people would want to test what they have down actually looks good and works well (especially links etc.) before pushing it to github. Luckily it is possible very easily: the docs website is powered by [mkdocs](https://www.mkdocs.org) which is a static website generator from markdown files (no database) with a yaml configuration file (look at the `mkdocs.yml` in the root of iLOVECLIM if interested) with a python framework on top of it. If you run a recent version of python, you should consider creating a virtual environnement and include the `mkdocs` package in it (something along the lines of `pip install mkdocs`). Once it is done, go to the iLOVECLIM root and run:

```bash
mkdocs serve
```
If you point a web-browser to the local address [http://127.0.0.1:8000](http://127.0.0.1:8000) you will see your local version of the website in action. If you keep it running like this, you will see that it updates live with your changes in the markdown files. Pretty cool, uh? :sunglasses:
