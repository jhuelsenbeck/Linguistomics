LIBS     =  -stdlib=libc++

CFLAGS   =  -Ofast -std=c++17

CC       =   clang++

OBJECTS  =  main.o Alignment.o AlignmentDistribution.o McmcSummary.o Msg.o Node.o ParameterStatistics.o Partition.o RandomVariable.o RbBitSet.o Sequence.o Subset.o Tree.o UserSettings.o Container.o

PROGS    = ttread

all:		$(PROGS)

ttread:		$(OBJECTS)
		$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o ttread
		
main.o:	main.cpp
		$(CC) $(CFLAGS) -c main.cpp

Alignment.o:	Alignment.cpp
		$(CC) $(CFLAGS) -c Alignment.cpp

AlignmentDistribution.o:	AlignmentDistribution.cpp
		$(CC) $(CFLAGS) -c AlignmentDistribution.cpp
		
Container.o:	Container.cpp
		$(CC) $(CFLAGS) -c Container.cpp
		
McmcSummary.o:	McmcSummary.cpp
		$(CC) $(CFLAGS) -c McmcSummary.cpp
		
Msg.o:	Msg.cpp
		$(CC) $(CFLAGS) -c Msg.cpp

Node.o:	Node.cpp
		$(CC) $(CFLAGS) -c Node.cpp

ParameterStatistics.o:	ParameterStatistics.cpp
		$(CC) $(CFLAGS) -c ParameterStatistics.cpp

Partition.o:	Partition.cpp
		$(CC) $(CFLAGS) -c Partition.cpp

RandomVariable.o:	RandomVariable.cpp
		$(CC) $(CFLAGS) -c RandomVariable.cpp

RbBitSet.o:	RbBitSet.cpp
		$(CC) $(CFLAGS) -c RbBitSet.cpp

Sequence.o:	Sequence.cpp
		$(CC) $(CFLAGS) -c Sequence.cpp

Subset.o:	Subset.cpp
		$(CC) $(CFLAGS) -c Subset.cpp

Tree.o:	Tree.cpp
		$(CC) $(CFLAGS) -c Tree.cpp

UserSettings.o:	UserSettings.cpp
		$(CC) $(CFLAGS) -c UserSettings.cpp

clean:		
		rm -f *.o
