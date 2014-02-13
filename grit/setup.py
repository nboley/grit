from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

config = {
    'cmdclass': {'build_ext': build_ext},
    'ext_modules': [
        Extension("sparsify_support_fns", ["grit/sparsify_support_fns.pyx"]),
    ],    
    'description': 'GRIT',
    'author': 'Nathan Boley',
    'url': 'http://grit-bio.org/',
    'download_url': 'http://grit-bio.org/git/',
    'author_email': 'npboley@gmail.com',
    'version': '1.0.0',
    'packages': ['grit', 'grit.analyze', 'grit.build', 'grit.files', 
                 'grit.lib', 'grit.proteomics', 
                 'grit.random_forest', 'grit.simulator'],
                 
    'setup_requires': ['Cython','numpy','scipy', 'networkx', 'igraph', 'pysam'],
    'scripts': [],
    'name': 'GRIT'
}

for_ubuntu = """
sudo aptitude install cython  python-dev python-numpy python-scipy \
                      python-igraph python-pip python-networkx
easy_install pysam
"""

def verify_dependencies():
    try: import scipy
    except: raise ValueError, "GRIT requires scipy"
    try: import igraph
    except: raise ValueError, "GRIT requires igraph"
    try: import networkx
    except: raise ValueError, "GRIT requires networkx"
    try: import pysam
    except: raise ValueError, "GRIT requires pysam"
    

if __name__== '__main__':
    verify_dependencies()
    setup(**config)
