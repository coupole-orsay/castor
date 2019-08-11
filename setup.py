import setuptools

with open('README.md', 'r') as f:
    long_description = f.read()
with open('requirements.txt', 'r') as f:
    requirements = f.read().strip('\n').split('\n')

entry_points = {
    'console_scripts': [
        'castor_exoplanet_analysis = castor._console_tools.exoplanet_analysis:main',
        'castor_pointing_analysis = castor._console_tools.pointing_analysis:main',
        'castor_prepare = castor._console_tools.prepare:main',
        'castor_align = castor._console_tools.align:main',
        ]
    }

package_data = {
    '': ['data/*'],
    }

setuptools.setup(
    name='castor',
    version='2019.08.08',
    author='Gabriel Pelouze',
    author_email='gabriel.pelouze@ias.u-psud.fr',
    description='Codes pour l’ASTronomie à ORsay',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/gpelouze/castor',
    packages=setuptools.find_packages(),
    entry_points=entry_points,
    package_data=package_data,
    python_requires='>=3.5',
    install_requires=requirements,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
)
