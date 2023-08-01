OBJDIR=../obj
BINDIR=../bin
LIBDIR=../lib

include ../Makerules

mulAcst: $(BINDIR)/mulAcst
	
$(BINDIR)/mulAcst:   \
	$(OBJDIR)/mulAcst.o
	$(FLINK) -o $@ \
	$(OBJDIR)/mulAcst.o \

$(OBJDIR)/mulAcst.o: mulAcst.f
	$(FCOMPILE) -o $@ mulAcst.f



misfit-permode2: $(BINDIR)/misfit-permode2
	
$(BINDIR)/misfit-permode2:   \
	$(OBJDIR)/misfit-permode2.o
	$(FLINK) -o $@ \
	$(OBJDIR)/misfit-permode2.o \

$(OBJDIR)/misfit-permode2.o: misfit-permode2.f
	$(FCOMPILE) -o $@ misfit-permode2.f



avmisfit-permode2: $(BINDIR)/avmisfit-permode2
	
$(BINDIR)/avmisfit-permode2:   \
	$(OBJDIR)/avmisfit-permode2.o
	$(FLINK) -o $@ \
	$(OBJDIR)/avmisfit-permode2.o \

$(OBJDIR)/avmisfit-permode2.o: avmisfit-permode2.f
	$(FCOMPILE) -o $@ avmisfit-permode2.f



compu-deriv-dstmodel-bspl: $(BINDIR)/compu-deriv-dstmodel-bspl
	
$(BINDIR)/compu-deriv-dstmodel-bspl:   \
	$(OBJDIR)/compu-deriv-dstmodel-bspl.o 
	$(FLINK) -o $@ \
	$(OBJDIR)/compu-deriv-dstmodel-bspl.o \
	$(LIBDIR)/libandrew.a \
	$(LIBDIR)/libgeo.a 

$(OBJDIR)/compu-deriv-dstmodel-bspl.o: compu-deriv-dstmodel-bspl.f
	$(FCOMPILE) -o $@ compu-deriv-dstmodel-bspl.f



m2sph-6l: $(BINDIR)/m2sph-6l
	
$(BINDIR)/m2sph-6l:   \
	$(OBJDIR)/m2sph-6l.o
	$(FLINK) -o $@ \
	$(OBJDIR)/m2sph-6l.o  

$(OBJDIR)/m2sph-6l.o: m2sph-6l.f
	$(FCOMPILE) -o $@ m2sph-6l.f




sph2m-6l: $(BINDIR)/sph2m-6l
	
$(BINDIR)/sph2m-6l:   \
	$(OBJDIR)/sph2m-6l.o
	$(FLINK) -o $@ \
	$(OBJDIR)/sph2m-6l.o  

$(OBJDIR)/sph2m-6l.o: sph2m-6l.f
	$(FCOMPILE) -o $@ sph2m-6l.f



bspl2v-3l: $(BINDIR)/bspl2v-3l
	
$(BINDIR)/bspl2v-3l:   \
	$(OBJDIR)/bspl2v-3l.o
	$(FLINK) -o $@ \
	$(OBJDIR)/bspl2v-3l.o  

$(OBJDIR)/bspl2v-3l.o: bspl2v-3l.f
	$(FCOMPILE) -o $@ bspl2v-3l.f



invATA-lapack-layer: $(BINDIR)/invATA-lapack-layer
	
$(BINDIR)/invATA-lapack-layer:   \
	$(OBJDIR)/invATA-lapack-layer.o 
	$(FLINK) -o $@ \
	$(OBJDIR)/invATA-lapack-layer.o \
	$(LIBDIR)/liblapack.a \
	$(LIBDIR)/libblas.a 

$(OBJDIR)/invATA-lapack-layer.o: invATA-lapack-layer.f
	$(FCOMPILE) -o $@ invATA-lapack-layer.f



buildATAcst: $(BINDIR)/buildATAcst
	
$(BINDIR)/buildATAcst:   \
	$(OBJDIR)/buildATAcst.o 
	$(FLINK) -o $@ \
	$(OBJDIR)/buildATAcst.o \
	$(LIBDIR)/libandrew.a \
	$(LIBDIR)/libgeo.a 

$(OBJDIR)/buildATAcst.o: buildATAcst.f
	$(FCOMPILE) -o $@ buildATAcst.f




addATAm: $(BINDIR)/addATAm
	
$(BINDIR)/addATAm:   \
	$(OBJDIR)/addATAm.o
	$(FLINK) -o $@ \
	$(OBJDIR)/addATAm.o \

$(OBJDIR)/addATAm.o: addATAm.f
	$(FCOMPILE) -o $@ addATAm.f
