#!/bin/bash
# If it has Build.PL use that, otherwise use Makefile.PL
if [ -f Build.PL ]; then
    perl Build.PL
    sed -i.bak -e '1 s|^.*$|#!/usr/bin/env perl|' Build
    ./Build installdeps
    # Make sure this goes in site
    ./Build install --installdirs site
elif [ -f Makefile.PL ]; then
    # Make sure this goes in site
    perl Makefile.PL INSTALLDIRS=site
    make
    make install
else
    echo 'Unable to find Build.PL or Makefile.PL. You need to modify build.sh.'
    exit 1
fi

export PERL5LIB="/opt/conda/lib/perl5/site_perl/5.22.0"
export PERL5LIB="/opt/conda/lib/perl5/site_perl/5.22.0/Modules"
export PERL5LIB="/opt/conda/lib/perl5/site_perl/5.22.0/x86_64-linux-thread-multi:$PERL5LIB"
# Add more build steps here, if they are necessary.

# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
