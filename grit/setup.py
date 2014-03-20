#from distutils.core import setup
#from distutils.extension import Extension

from setuptools import setup, Extension

try: import Cython
except ImportError: USE_CYTHON = False
else: USE_CYTHON = True

# cython seems to be broken on some architecture, so just use
# the c file until I can resolve this
USE_CYTHON = False
ext = 'pyx' if USE_CYTHON else 'c'

extensions = [
    Extension("sparsify_support_fns", ["grit/sparsify_support_fns." + ext, ])
]
if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

config = {
    'ext_modules': extensions,
    'description': 'GRIT',
    'author': 'Nathan Boley',
    'url': 'http://grit-bio.org/',
    'download_url': 'http://grit-bio.org/git/',
    'author_email': 'npboley@gmail.com',
    'version': '1.1.0',
    'packages': ['grit', 
                 'grit.analyze', 
                 'grit.files', 
                 'grit.lib', 
                 'grit.proteomics', 
                 'grit.random_forest', 
                 'grit.simulator'],
                 
    'setup_requires': [ 'scipy', 'numpy', 'Cython', 'networkx', 'pysam' ],
    'install_requires': [ 'scipy', 'numpy', 'Cython', 'networkx', 'pysam' ],
    'scripts': ['./bin/run_grit.py',],
    'name': 'GRIT'
}

if __name__== '__main__':
    setup(**config)
