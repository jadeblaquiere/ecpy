from setuptools import setup
from ecpy import __version__

setup(name='ECPY: Pure python elliptic curve math', 
      version=str(__version__),
      author='Joseph deBlaquiere',
      author_email='jadeblaquiere@yahoo.com',
      packages=['ecpy'])
