from setuptools import setup, find_packages
import os

# Read file content
def read_file(filename):
    filepath = os.path.join(os.path.dirname(__file__), filename)
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            return f.read()
    except FileNotFoundError:
        return ""

# Get version number
def get_version():
    init_path = os.path.join('qcombo', '__init__.py')
    if os.path.exists(init_path):
        with open(init_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('__version__'):
                    return line.split('=')[1].strip().strip("'\"")
    return '0.1.0'

# Read README and LICENSE
long_description = read_file('README.md')
license_text = read_file('LICENSE')

setup(
    # Basic information
    name='qcombo',
    version=get_version(),
    author='L.H.chen, Y.Li, Heiko Hergert, J.M.Yao',
    author_email='2942445940@qq.com',
    description='Quantum-Commutator of Many-Body Operators based on generalized Wick theorem',
    
    # Detailed description
    long_description=long_description,
    long_description_content_type='text/markdown',
    
    # License
    license='MIT',
    license_files=['LICENSE'],
    
    # Project URLs
    url='https://github.com/chenlh73/qcombo',
    
    # Package management
    packages=find_packages(exclude=['tests*', 'examples*', 'results*']),
    include_package_data=True,
    
    # Dependencies
    python_requires='>=3.12',
    install_requires=[
        'sympy>=1.13.3',
    ],

    entry_points={
        'console_scripts': [
            'qcombo=qcombo.__main__:cli_main',
        ],
    },

    # Keywords
    keywords=[
        'quantum',
        'wick-theorem',
        'commutator',
    ],

    classifiers=[
        # Development status
        'Development Status :: 3 - Alpha',
        
        # Target audience
        'Intended Audience :: Science/Research',
        
        # Topics
        'Topic :: Scientific/Engineering',
        
        # License
        'License :: OSI Approved :: MIT License',
        
        # Python versions (based on actual testing)
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.12',
        
        # Operating system
        'Operating System :: OS Independent',
    ],
)

