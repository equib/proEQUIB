from distutils.core import setup
from distutils.command.install import install
import json
import os.path
import sys

kernel_json = {"argv":[sys.executable,"-m","gdl_kernel","-f","{connection_file}"],
               "display_name":"GDL",
               "language":"GDL",
               "codemirror_mode":"gdl"}

class install_with_kernelspec(install):
    def run(self):
        install.run(self)

        from jupyter_client.kernelspec import KernelSpecManager
        from IPython.utils.path import ensure_dir_exists
        destdir = os.path.join(KernelSpecManager().user_kernel_dir,'GDL')
        ensure_dir_exists(destdir)
        with open(os.path.join(destdir,'kernel.json'),'w') as f:
            json.dump(kernel_json,f,sort_keys=True)


with open('gdl_kernel_readme.md') as f:
    readme = f.read()

svem_flag = '--single-version-externally-managed'
if svem_flag in sys.argv:
    sys.argv.remove(svem_flag)

setup(name='GDL_kernel',
      version='0.5',
      description='A GDL kernel for Jupyter',
      long_description=readme,
      author='Luke Stagner, GDL Team',
      url='https://github.com/gnudatalanguage/gdl_kernel',
      py_modules=['gdl_kernel'],
      cmdclass={'install': install_with_kernelspec},
#      install_requires=['pexpect>=3.3','IPython >= 3.0'],
      classifiers = [
          'Framework :: Jupyter',
          'License :: OSI Approved :: BSD License',
          'Programming Language :: Python :: 3',
          'Topic :: System :: Shells',
      ]
)
