######################################################################
# Automatically generated by qmake (1.07a) Wed Jan 13 16:14:23 2010
######################################################################

TEMPLATE = lib
CONFIG -= moc qt
DEPENDPATH += ./src
INCLUDEPATH += . ./src

unix{
   INCLUDEPATH += /home/pompeu/local64/include
   INCLUDEPATH += /usr/include
   LIBS += -L/home/pompeu/local64/lib
}
macx{
   CONFIG -= app_bundle
   INCLUDEPATH += /Users/mtcheou/local/include
   LIBS += -L/Users/mtcheou/local/lib
}

QMAKE_CFLAGS_RELEASE   =  -w -pipe -D_FORTIFY_SOURCE=2 \
                         -fomit-frame-pointer
QMAKE_CXXFLAGS_RELEASE =  -w -pipe -D_FORTIFY_SOURCE=2 \
                         -fomit-frame-pointer

# Input
HEADERS += ./src/bincomtd.h \
           ./src/buftipos.h \
           ./src/comtrade.h \
           ./src/newcomtd.h \
           ./src/newcomtd_99.h 
SOURCES += ./src/bincomtd.cpp \
           ./src/buftipos.cpp \
           ./src/comtrade.cpp \
           ./src/newcomtd.cpp \
           ./src/newcomtd_99.cpp

OBJECTS_DIR  = obj
#DESTDIR = ../lib
TARGET       = comtrade

header_files.files = $$HEADERS
header_files.path = ../include
INSTALLS += header_files

target.path = ../lib
INSTALLS += target
 
