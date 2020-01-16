#-------------------------------------------------
#
# Project created by QtCreator 2015-04-24T09:42:56
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = MBL2017
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    lineage.cpp \
    simulation.cpp \
    genus.cpp \
    qcustomplot.cpp \
    about.cpp \
    darkstyletheme.cpp

HEADERS  += mainwindow.h \
    lineage.h \
    simulation.h \
    genus.h \
    qcustomplot.h \
    version.h \
    about.h \
    darkstyletheme.h

FORMS    += mainwindow.ui \
    about.ui

RESOURCES += resources.qrc

#Mac icon
ICON = ./resources/mbl.icns

#Needed to make binaries launchable from file in Ubuntu - GCC default link flag -pie on newer Ubuntu versions this so otherwise recognised as shared library
QMAKE_LFLAGS += -no-pie
