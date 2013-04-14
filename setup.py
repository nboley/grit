try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'GRIT',
    'author': 'Nathan Boley',
    'url': 'http://grit-bio.org/',
    'download_url': 'http://grit-bio.org/git/',
    'author_email': 'npboley@gmail.com',
    'version': '1.0',
    'install_requires': [''],
    'packages': ['grit'],
    'scripts': [],
    'name': 'GRIT'
}

setup(**config)

