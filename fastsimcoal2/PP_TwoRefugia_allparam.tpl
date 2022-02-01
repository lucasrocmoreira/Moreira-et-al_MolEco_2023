//Parameters for the coalescence simulation program : fastsimcoal.exe
4 samples to simulate :
//Population effective sizes (number of genes)
NPOP1
NPOP2
NPOP3
NPOP4
//Samples sizes and samples age (0=AK, 1=E, 2=R, 3=NW)
10 
10
10
10
//Growth rates	: negative growth implies population expansion
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0	M21	M31	M41
M12	0	M32	M42
M13	M23	0	M43
M14	M24	M34	0
//Migration matrix 1
0	0	0	0
0	0	0	0
0	0	0	0
0	0	0	0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
7 historical event
tch1 0 0 1 SC1 0 1
tch2 1 1 1 SC2 0 1
tch3 2 2 1 SC3 0 1
tch4 3 3 1 SC4 0 1
tchDiff1 3 1 1 RESIZE1 0 1
tchDiff2 2 1 1 RESIZE2 0 1
tchDiff3 0 1 1 RESIZE3 0 1
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0.0000 2.42e-9 OUTEXP



























