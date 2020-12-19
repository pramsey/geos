#
# This file is part of project GEOS (http://trac.osgeo.org/geos/)
#
prefix=@prefix@
top_srcdir=@top_srcdir@
top_builddir=@top_builddir@

DOXYGEN = doxygen

doxygen-html apidoc doxygen: Doxyfile
	$(DOXYGEN) $<

clean-local:
	rm -Rf doxygen_docs

EXTRA_DIST = Doxyfile \
	CMakeLists.txt
