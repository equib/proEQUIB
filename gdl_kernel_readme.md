# GDL kernel for Jupyter

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/equib/gdl_kernel/HEAD?labpath=notebooks%2Fdemo_gdl.ipynb)

Demo: [Notebook](https://github.com/equib/gdl_kernel/blob/master/notebooks/demo_gdl.ipynb) on [Binder](https://mybinder.org/v2/gh/equib/gdl_kernel/HEAD?labpath=notebooks%2Fdemo_gdl.ipynb).

This package is originally from [github.com/gnudatalanguage/gdl_kernel](https://github.com/gnudatalanguage/gdl_kernel) and modified to work with [Binder](https://mybinder.org/). This repository also makes use of [GDL](https://github.com/gnudatalanguage/gdl) to run IDL/GDL codes in Jupyter lab. This requires the following component:

* `apt.txt` for apt-installing the gnudatalanguage (gdl) component.

The current version was found to work with pexpect 4.6, jupyter 4.4, and jupyter-notebook 5.7. 

To install:
```
python setup.py install --prefix=~/.local
```
or
```
python3 setup.py install --prefix=~/.local
```

This should make an GDL directory containing the kernelspec in `~/.local/share/jupyter/kernels`, and place `gdl_kernel.py` in `~/.local/lib/pythonX.X/site-packages`.

To run:
```
jupyter notebook 
#Select GDL kernel from dropdown menu
```
