from setuptools import setup, Extension, find_packages

extensions = [
    Extension("grit.sparsify_support_fns", ["grit/sparsify_support_fns.pyx", ]),
    Extension("grit.call_peaks_support_fns", ["grit/call_peaks_support_fns.pyx", ])\
]

config = {
    'include_package_data': True,
    'ext_modules': extensions, 
    'description': 'GRIT',
    'author': 'Nathan Boley',
    'url': 'http://grit-bio.org/',
    'download_url': 'http://grit-bio.org/git/',
    'author_email': 'npboley@gmail.com',
    'version': '2.0.0',
    'packages': ['grit', 
                 'grit.analyze', 
                 'grit.files', 
                 'grit.lib', 
                 'grit.proteomics'],
    'setup_requires': [ 'Cython', ],
    'install_requires': [ 'scipy', 'numpy', 'Cython', 'networkx', 'pysam' ],
    'scripts': ['./bin/run_grit.py', "./bin/call_peaks.py"],
    'name': 'GRIT'
}

if __name__== '__main__':
    setup(**config)
