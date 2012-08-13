#this is the makefile for GC class
RM=rm -f
TAR=tar
JARLIBHOME=./../lib
CLASSPATH=$(JARLIBHOME)/jar/*
JAVA=java
#DEBUGFLAGS= -g 
JAVAFLAGS=-cp $(CLASSPATH):bin $(DEBUGFLAGS)
JAVAC=javac -d bin -sourcepath src
JAVACFLAGS=-cp $(CLASSPATH):. -encoding utf-8 -Xlint:deprecation

SOURCES=src/SimpleExampleBase.java src/SimpleExample.java src/Icosahedron.java src/Object3D.java src/ImageLoader.java src/Shader.java src/DotImage.java src/TextureImage.java src/TextureImageUtil.java
#($SOURCES:%=bin/%)

.PHONY: clean cleanall remote tar
.SUFFIXES : .java .class

ALL::   bin/SimpleExample.class
	$(JAVA) $(JAVAFLAGS) $(<:bin/%.class=%) 

bin/SimpleExample.class: $(SOURCES)
	mkdir -p bin
	$(JAVAC) $(JAVACFLAGS) $(@:bin/%.class=src/%.java)

clean:
	$(RM) -R bin

cleanall:
	$(RM) -R bin javadoc

tar:
	$(RM) -R bin javadoc
	cd ..; $(TAR) zcvf simpleexample.tar.gz SimpleExample

remote: bin/SimpleExample.class
#	DISPLAY=localhost:0.0;$(JAVA) $(<:bin/%.class=%) &
	DISPLAY=:0.0;$(JAVA) $(JAVAFLAGS) $(<:bin/%.class=%) &
	sleep 5
#       DISPLAY=localhost:0.0;xwd -name simpleJoglShader -out t
#       DISPLAY=localhost:0.0;xwd -root -out t
	DISPLAY=:0.0;xwd -name 3DCGSimpleExample -out t
	xwud < t
	rm -f t
	killall java
