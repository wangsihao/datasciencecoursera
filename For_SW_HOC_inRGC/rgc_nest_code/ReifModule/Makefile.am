# Automake file for external dynamic modules for NEST
#
# Hans Ekkehard Plesser, April 2008
# Automake file for the Developer Module
#
# libmymodule is built as a normal, installable library.
# It will be installed to $prefix/lib by make install.
#
# Headers from this directory are not to be installed upon
# make install. They are therefore included in _SOURCES.


# 1. Exchange "my" in "mymodule" with the name of your model below
#    (ten locations).
# 2. Add all .cpp and .h files from your code as *_SOURCES. Header files
#    are given only so that they will be included in the tarball if you
#    run "make dist" on your module.
# 3. The mymodule* stuff creates a module that can be loaded at runtime.
#    It is called mymodule.so.
# 4. The libmymodule* stuff creates a library against which NEST can be
#    linked.

libdir= @libdir@/nest

lib_LTLIBRARIES=      reifmodule.la libreifmodule.la

reifmodule_la_CXXFLAGS= @AM_CXXFLAGS@
reifmodule_la_SOURCES=  reifmodule.cpp      reifmodule.h     \
			reif_cond_alpha_RK5.cpp   reif_cond_alpha_RK5.h  \
                      pif_psc_alpha.cpp pif_psc_alpha.h \
                      drop_odd_spike_connection.h
reifmodule_la_LDFLAGS=  -module

libreifmodule_la_CXXFLAGS= $(reifmodule_la_CXXFLAGS) -DLINKED_MODULE
libreifmodule_la_SOURCES=  $(reifmodule_la_SOURCES)

MAKEFLAGS= @MAKE_FLAGS@

AM_CPPFLAGS= @NEST_CPPFLAGS@ \
             @INCLTDL@

.PHONY: install-slidoc

pkgdatadir=@datadir@/nest

nobase_pkgdata_DATA=\
	sli/reifmodule-init.sli

install-slidoc:
	NESTRCFILENAME=/dev/null $(DESTDIR)$(NEST_PREFIX)/bin/sli --userargs="@HELPDIRS@" $(NEST_PREFIX)/share/nest/sli/install-help.sli

install-data-hook: install-exec install-slidoc

EXTRA_DIST= sli
