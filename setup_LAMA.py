#!/usr/bin/env python

print "=== LAMA phenotyping pipeline ===\nDownloading dependencies..."


# easy_install first
try:
    from setuptools.command import easy_install
    dependencies = ["numpy", "scipy" "pyyaml", "SimpleITK"]
    dependencies_named = {
        "yaml": "pyyaml",

    }

    for dep in dependencies:

        try:
            print "Installing {0}...".format(dep),
            mod = __import__(dep)  # try to import module
            print " already installed.".format(dep)

        except ImportError:
            # If it fails, try to easy install it
            easy_install.main(["--user", dep])
            print "done."

except ImportError:
    print "Couldn't locate 'easy_install'. Do you have setuptools installed on your machine? Try sudo apt-get install python-setuptools (Linux) or use Homebrew (Mac)."

