#~/bin/sh
# Re-generate autoTools files
mkdir -p autoconf/
aclocal 
autoconf
autoheader
libtoolize
automake --add-missing

