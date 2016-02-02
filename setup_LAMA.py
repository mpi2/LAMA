#!/usr/bin/env python

print "=== LAMA phenotyping pipeline ===\nDownloading dependencies..."


# easy_install first
try:
    from setuptools.command import easy_install
    dependencies = ["numpy", "scipy", "SimpleITK", 'appdirs', 'psutil']
    dependencies = {
        'yaml': 'pyyaml',
        'scipy': 'scipy',
        'numpy': 'numpy',
        'SimpleITK': 'SimpleITK',
        'appdirs': 'appdirs',
        'psutil': 'psutil',
        'yaml': 'pyyaml'

    }

    for import_name, package_name in dependencies.iteritems():

        try:
            print "Installing {0}...".format(import_name),
            mod = __import__(import_name)  # try to import module
            print " already installed.".format(import_name)

        except ImportError:
            # If it fails, try to easy install it
            easy_install.main(["--user", package_name])
            print "done."

except ImportError:
    print "Couldn't locate 'easy_install'. Do you have setuptools installed on your machine? Try sudo apt-get install python-setuptools (Linux) or use Homebrew (Mac)."

