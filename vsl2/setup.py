# Copyright Christopher Oldfield

from distutils.core import setup

setup(
    name='vsl2',
    packages=[
        'vsl2',
        'vsl2.data'
    ],
    package_data={
        'vsl2.data' : ['VSL2.jar']
    }
)
