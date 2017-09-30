#-----------------------------------------------------------------------------
# Top-level Makefile for ANN.
#
# ANN: Approximate Nearest Neighbors
# Version: 1.1 05/03/05
#
# THIS MAKEFILE MODIFIED FROM THE ORIGINAL BY GHAZALE HARATINEZHAD
# TO ALLOW COMPILATION AS AN R PACKAGE.
#-----------------------------------------------------------------------------
# Copyright (c) 1997-2005 University of Maryland and Sunil Arya and
# David Mount.  All Rights Reserved.
# 
# This software and related documentation is part of the Approximate
# Nearest Neighbor Library (ANN).  This software is provided under
# the provisions of the Lesser GNU Public License (LGPL).  See the
# file ../ReadMe.txt for further information.
# 
# The University of Maryland (U.M.) and the authors make no
# representations about the suitability or fitness of this software for
# any purpose.  It is provided "as is" without express or implied
# warranty.
#-----------------------------------------------------------------------------
# Revision 0.1  09/06/97
#	alpha release
# Revision 0.2  06/26/98
#	Minor changes to fix compilation errors on SGI systems.
# Revision 1.0  04/01/05
#	Initial release (finally!)
#	Added linux-g++ target
# Revision 1.1  05/03/05
#	Added macosx-g++ target
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# main make entry point
#-----------------------------------------------------------------------------
all:
	cd src ; $(MAKE) 
	cd test ; $(MAKE) 
	cd sample ; $(MAKE) 
	cd ann2fig ; $(MAKE) 

#-----------------------------------------------------------------------------
# Remove .o files and core files
#-----------------------------------------------------------------------------
clean:
	cd src ; $(MAKE) clean
	cd test ; $(MAKE) clean
	cd sample ; $(MAKE) clean
	cd ann2fig ; $(MAKE) clean
	cd doc ; $(MAKE) clean
	cd MS_Win32; $(MAKE) clean
	cd validate ; $(MAKE) clean

#-----------------------------------------------------------------------------
# Remove everthing that can be remade
#-----------------------------------------------------------------------------
realclean: 
	-rm -f lib/*
	-rm -f bin/*
	cd src ; $(MAKE) realclean
	cd test ; $(MAKE) realclean
	cd sample ; $(MAKE) realclean
	cd ann2fig ; $(MAKE) realclean
	cd doc ; $(MAKE) realclean
	cd MS_Win32; $(MAKE) realclean
	cd validate ; $(MAKE) realclean


