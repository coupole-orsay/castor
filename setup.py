import setuptools

with open('README.md', 'r') as f:
    long_description = f.read()
with open('requirements.txt', 'r') as f:
    requirements = f.read().strip('\n').split('\n')

entry_points = {
    'console_scripts': [
        'castor_align = castor._console_scripts.align:main',
        'castor_exoplanet_analysis = castor._console_scripts.exoplanet_analysis:main',
        'castor_pointing_analysis = castor._console_scripts.pointing_analysis:main',
        'castor_prepare = castor._console_scripts.prepare:main',
        'castor_rotate_spectra = castor._console_scripts.rotate_spectra:main',
        'castor_wavelength_calibration = castor._console_scripts.wavelength_calibration:main',
        'castor_align_spectra = castor._console_scripts.align_spectra:main',
        ]
    }

package_data = {
    '': ['data/*'],
    }

setuptools.setup(
    name='castor-orsay',
    version='2019.12.25.1',
    author='Gabriel Pelouze, Aurélien Stcherbinine',
    author_email='gabriel.pelouze@ias.u-psud.fr, aurelien.stcherbinine@ias.u-psud.fr',
    description='Codes pour l’ASTronomie à ORsay',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/coupole-orsay/castor',
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
