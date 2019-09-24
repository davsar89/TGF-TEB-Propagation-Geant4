#-------------------------------------------------------------------------------
# /* GEANT4 code for propagation of gamma-rays, electron and positrons in Earth's atmosphere */
# 
# //
# // ********************************************************************
# // * License and Disclaimer                                           *
# // *                                                                  *
# // * The  Geant4 software  is  copyright of the Copyright Holders  of *
# // * the Geant4 Collaboration.  It is provided  under  the terms  and *
# // * conditions of the Geant4 Software License,  included in the file *
# // * LICENSE and available at  http://cern.ch/geant4/license .  These *
# // * include a list of copyright holders.                             *
# // *                                                                  *
# // * Neither the authors of this software system, nor their employing *
# // * institutes,nor the agencies providing financial support for this *
# // * work  make  any representation or  warranty, express or implied, *
# // * regarding  this  software system or assume any liability for its *
# // * use.  Please see the license in the file  LICENSE  and URL above *
# // * for the full disclaimer and the limitation of liability.         *
# // *                                                                  *
# // * This  code  implementation is the result of  the  scientific and *
# // * technical work of the GEANT4 collaboration.                      *
# // * By using,  copying,  modifying or  distributing the software (or *
# // * any work based  on the software)  you  agree  to acknowledge its *
# // * use  in  resulting  scientific  publications,  and indicate your *
# // * acceptance of all terms of the Geant4 Software license.          *
# // ********************************************************************
#-------------------------------------------------------------------------------
# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.9 2010-01-11 16:20:08 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := XrayTel
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

ifdef G4ANALYSIS_USE_PLOTTER
  CPPFLAGS += -DG4ANALYSIS_USE_PLOTTER
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

