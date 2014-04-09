from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules = [Extension("hello",["hello.pyx"])]

setup(
     name = 'Hello World app',
     cmdclass = {'build_ext':build_ext},
     ext_modules = [Extension("countSeq2", ["/hpc/users/ummata01/gitrepos/workIn/debug/pacBioTR/psCount/countSeq2.pyx"], include_dirs=[numpy.get_include()])]
)

