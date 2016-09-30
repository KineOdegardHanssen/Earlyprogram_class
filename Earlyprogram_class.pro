TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    systems.cpp \
    diagonalization.cpp

HEADERS += \
    systems.h \
    diagonalization.h \

LIBS += -larmadillo -llapack -lblas
QMAKE_CXXFLAGS += -Wall -std=c++0x
INCLUDEPATH += "/home/ubu/Downloads/eigen-eigen-dc6cfdf9bcec/"
INCLUDEPATH += "/home/ubu/Downloads/lapack-3.6.1/LAPACKE/include/"
INCLUDEPATH += "/home/ubu/Downloads/lapack-3.6.1/"

