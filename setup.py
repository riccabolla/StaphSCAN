from setuptools import setup, find_packages

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="staphscan",
    version="0.3.0",  
    description="A tool for Staphylococcus aureus analysis",
    author="Riccardo Bollini",  
    url="https://github.com/riccabolla/StaphSCAN", 
    
    packages=find_packages(),
        
    package_data={
        'staphscan': ['modules/*/data/*']
    },
    include_package_data=True,
    
    install_requires=requirements,
    
    entry_points={
        'console_scripts': [
            'staphscan=staphscan.main:main',  
        ],
    },
    
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
)