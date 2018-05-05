#!/usr/bin/env python

print "=== LAMA phenotyping pipeline ===\nDownloading dependencies..."


# easy_install first
import sys
try:
    import pip
except ImportError:
    print "setup_LAMA requires 'pip' to be installed on your system\non ubuntu try 'sudo apt install python-pip'"
    sys.exit()

dependencies = {
    'scipy': 'scipy',
    'numpy': 'numpy',
    'SimpleITK': 'SimpleITK',
    'appdirs': 'appdirs',
    'psutil': 'psutil',
    'yaml': 'pyyaml',
    'sklearn': 'sklearn',
    'matplotlib': 'matplotlib',
    'pandas': 'pandas',
    'seaborn': 'seaborn',
    'pandas': 'pandas',
    'statsmodels': 'statsmodels',
    'PIL': 'Pillow'

}

failed_installs = []

for import_name, package_name in dependencies.iteritems():

    try:
        print "Installing {0}...".format(import_name),
        mod = __import__(import_name)  # try to import module
        print " already installed.".format(import_name)

    except ImportError:
        # If it fails, try to install with pip
        result = pip.main(['install', '--user', package_name])
        if result != 0:
            failed_installs.append(package_name)

if len(failed_installs) == 0:
    print "All packages successfully installed"
else:
    print 'The following packages failed to install\n'
    for failed in failed_installs:
        print "{}\n".format(failed)



