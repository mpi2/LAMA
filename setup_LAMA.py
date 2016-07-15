#!/usr/bin/env python

print "=== LAMA phenotyping pipeline ===\nDownloading dependencies..."


# easy_install first
try:
    import pip
except ImportError:
    print "setup_LAMA requires 'pip' to be installed on your system\n"

dependencies = {
    'yaml': 'pyyaml',
    'scipy': 'scipy',
    'numpy': 'numpy',
    'SimpleITK': 'SimpleITK',
    'appdirs': 'appdirs',
    'psutil': 'psutil',
    'yaml': 'pyyaml'

}

failed_installs = []

for import_name, package_name in dependencies.iteritems():

    try:
        print "Installing {0}...".format(import_name),
        mod = __import__(import_name)  # try to import module
        print " already installed.".format(import_name)

    except ImportError:
        # If it fails, try to easy install it
        result = pip.main(["--user", package_name])
        if result != 0:
            failed_installs.append(package_name)

if len(failed_installs) == 0:
    print "All packages successfully installed"
else:
    print 'The following packages failed to install\n'
    for failed in failed_installs:
        print "{}\n".format(failed)

if 'SimpleITK' in failed_installs:

    itk_url = 'http://sourceforge.net/projects/simpleitk/files/SimpleITK/0.9.0/Python/SimpleITK-0.9.0-py2.7-linux-x86_64.egg'
    print "\nSometimes SimpleITK can be problematic to install\n"
    print "Download the puthon egg from here\n\t{}\n\n".format(itk_url)
    print "And install with this command\n"
    print "easy_install --user SimpleITK-0.9.0-py2.7-linux-x86_64.egg"


