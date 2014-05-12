(* ::Package:: *)

DataCreate[distribution_,size_]:=
{
data=RandomVariate[distribution,size],
dataplot=ListPlot[data],
datapdfplot=Plot[PDF[NormalDistribution[],x],{x,-5,5}]
}


DataCreate[NormalDistribution[],1000];
datapdfplot


Kerneling[]:=
{
kerdist=SmoothKernelDistribution[data],
kerpdfplot=Plot[PDF[kerdist,x],{x,-5,5}],
kercdfplot=Plot[CDF[kerdist,x],{x,-5,5}]
}


Kerneling[];
kerpdfplot
kercdfplot



