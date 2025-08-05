SUBDIR = \
	Kurosh Beam Buffer Contact ContactFast \
	Libs Link Material Timing EB_Beam EB_Beam_A

all:	mumps sparse math submodules build

submodules: distrib
	for d in $(SUBDIR) ; \
		do \
			echo $$d xxxxxxxxxxxxxx ; \
			cd $$d ; ./ANSCUSTOM ; cd .. ; \
		done
mumps:
	cp -f ./ANSCUSTOM_MUMPS MUMPS
	cd MUMPS ; ./ANSCUSTOM_MUMPS

math:
	cp -f ./ANSCUSTOM_MATH Math
	cd Math ; ./ANSCUSTOM_MATH

sparse:
	cp -f ./ANSCUSTOM_SPARSE  Sparse
	cd Sparse ; ./ANSCUSTOM_SPARSE

build:
	rm -f libAnsEB_Beam.a
	rm -f libAnsEB_Beam_A.a
	./ANSCUSTOM

distrib:
	for d in $(SUBDIR) ; \
		do \
			cp -f ./ANSCUSTOM_NOLINK $$d/ANSCUSTOM ; \
		done
