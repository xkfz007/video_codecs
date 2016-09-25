TARGET:=attention
INCDIR:=inc
SRCDIR:=src
OBJDIR:=obj
$(shell if !([ -d $(OBJDIR) ]);then mkdir $(OBJDIR);fi)
LIBDIR:=lib
$(shell if !([ -d $(LIBDIR) ]);then mkdir $(LIBDIR);fi)

CC  :=g++

FLAGS:=-Wall -W -I$(INCDIR) -D_DEBUG
LIBS :=-lm

OPENCV=1 #switch of opencv
ifdef OPENCV
FLAGS+=`pkg-config --cflags opencv`
LIBS +=`pkg-config --libs opencv`
endif

#DBG=0
ifdef DBG
FLAGS+=-g
else
FLAGS+=-O2
endif

SRC:=$(wildcard $(SRCDIR)/*.cpp)
OBJ:=$(patsubst %.cpp,$(OBJDIR)/%.o,$(notdir $(SRC)))

#all:
#	@echo $(OBJ)
#

$(TARGET):$(OBJ)
	$(CC) -o $@ $^ $(FLAGS) $(LIBS)
$(OBJDIR)/%.o:$(SRCDIR)/%.cpp 
	$(CC) -c $< -o $@ $(FLAGS)
.PHONY: clean
clean:
	rm -f $(OBJDIR)/*.o *~ core 

