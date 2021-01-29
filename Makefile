CXX = g++
CXXFLAGS  += -O2 -Wall -Werror -Wextra
RCXXFLAGS := `root-config --cflags --libs`
LDFLAGS   += -lconf -L./git/config/lib \
	     -lfoliage -L./git/foliage/lib \
	     -lhist -L./git/history/lib \
	     -lpp -L./git/paper-and-pencil/lib \
	     -ltt -L./git/tricks-and-treats/lib
RLDFLAGS  := -lEG -lTMVA -lUnfold

BINDIR = ./bin
BLDDIR = ./build
SRCDIR = ./src
TSTDIR = ./test

RPPSRCS = $(wildcard $(SRCDIR)/*.C)
RPPEXES = $(patsubst $(SRCDIR)/%.C,$(BINDIR)/%,$(RPPSRCS))
RPPDEPS = $(patsubst $(SRCDIR)/%.C,$(BLDDIR)/%.d,$(RPPSRCS))

CPPSRCS = $(wildcard $(SRCDIR)/*.cpp)
CPPEXES = $(patsubst $(SRCDIR)/%.cpp,$(BINDIR)/%,$(CPPSRCS))
CPPDEPS = $(patsubst $(SRCDIR)/%.cpp,$(BLDDIR)/%.d,$(CPPSRCS))

TSTSRCS = $(wildcard $(TSTDIR)/*.C)
TSTEXES = $(patsubst $(TSTDIR)/%.C,$(BINDIR)/%,$(TSTSRCS))
TSTDEPS = $(patsubst $(TSTDIR)/%.C,$(BLDDIR)/%.d,$(TSTSRCS))

EXES = $(RPPEXES) $(CPPEXES) $(TSTEXES)
DEPS = $(RPPDEPS) $(CPPDEPS) $(TSTDEPS)

.PHONY: all submodules clean

all: $(RPPEXES) $(CPPEXES) $(TSTEXES)

submodules:
	$(MAKE) -C ./git/config
	$(MAKE) -C ./git/foliage
	$(MAKE) -C ./git/history
	$(MAKE) -C ./git/paper-and-pencil
	$(MAKE) -C ./git/tricks-and-treats

$(BINDIR)/%: $(SRCDIR)/%.C
	@mkdir -p $(BINDIR)
	@mkdir -p $(BLDDIR)
	$(CXX) $(CXXFLAGS) $(RCXXFLAGS) -MMD -MF $(BLDDIR)/$(*F).d $< -o $@ \
		$(LDFLAGS) $(RLDFLAGS)

$(BINDIR)/%: $(SRCDIR)/%.cpp
	@mkdir -p $(BINDIR)
	@mkdir -p $(BLDDIR)/$(@D)
	$(CXX) $(CXXFLAGS) -MMD -MF $(BLDDIR)/$(*F).d $< -o $@ \
		$(LDFLAGS)

$(BINDIR)/%: $(TSTDIR)/%.C
	@mkdir -p $(BINDIR)
	@mkdir -p $(BLDDIR)
	$(CXX) $(CXXFLAGS) $(RCXXFLAGS) -MMD -MF $(BLDDIR)/$(*F).d $< -o $@ \
		$(LDFLAGS) $(RLDFLAGS)

clean:
	@$(RM) $(EXES) $(DEPS)
	@rm -f $(BINDIR)/*
	@rm -rf $(BLDDIR)/*

-include $(DEPS)
