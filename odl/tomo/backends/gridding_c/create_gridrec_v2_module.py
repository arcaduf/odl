from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("gridrec_v2",
                             sources=[ "gridrec_v2.pyx" ,
                                       "gridrec_v2_backproj.c" ,
                                       "gridrec_v2_forwproj.c" ,
                                       "pswf.c" ,
                                       "filters.c"],
                             include_dirs=[numpy.get_include()],libraries=['fftw3f','gcov'],extra_compile_args=['-O3','-march=native','-ffast-math','-fprofile-generate'],extra_link_args=['-fprofile-generate'])],
)

'''
import gridrec_v2
gridrec_v2.createFFTWWisdomFile(2016, "profile.wis")

import os
import sys
os.system(sys.executable + " profile.py")

os.remove('gridrec_v2.so')

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("gridrec_v2",
                             sources=[ "gridrec_v2.pyx" ,
                                       "gridrec_v2_backproj.c" ,
                                       "gridrec_v2_forwproj.c" ,
                                       "pswf.c" ,
                                       "filters.c"],
                             include_dirs=[numpy.get_include()],libraries=['fftw3f'],extra_compile_args=['-O3','-march=native','-ffast-math','-fprofile-use'],extra_link_args=['-fprofile-use'])],
)
'''
