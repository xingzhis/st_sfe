from setuptools import setup, find_packages

setup(
    name='stfe',
    version='0.1.1',
    # packages=find_packages(where='stfe'),
    packages=find_packages(),
    # package_dir={'': 'stfe'},
    include_package_data=True,
    package_data={
        '': ['reference/*.xlsx', 'reference/*.txt'],
    },
    author="Xingzhi Sun",
    author_email="xingzhi.sun@yale.edu",
    description="Your package description",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    install_requires=[
        # List your package dependencies here
        # 'some_package>=1.0',
    ],
)
