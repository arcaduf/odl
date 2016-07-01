import os
import shutil


if os.path.exists('build'):
    shutil.rmtree('build')

os.system('python create_gridrec_v2_module.py build_ext --inplace')
os.rename( 'gridrec_v2.cpython-35m-x86_64-linux-gnu.so' , 'gridrec_v2.so' )
