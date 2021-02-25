#include "itensor/all.h"
#include "itensor/util/stdx.h"
using namespace itensor;
void calculatelocal(MPS const& psio,  SiteSet const& sites)
{    
	auto N=sites.N();
	Vector upd(N),dnd(N),szi(N);
	auto psi=psio;
	for(int j = 1; j <= N; ++j)
        {
        psi.position(j);
        upd(j-1) = (dag(prime(psi.A(j),Site))*sites.op("Nup",j)*psi.A(j)).real();
        dnd(j-1) = (dag(prime(psi.A(j),Site))*sites.op("Ndn",j)*psi.A(j)).real();
        szi(j-1) = (dag(prime(psi.A(j),Site))*sites.op("Sz",j)*psi.A(j)).real();
        }

    println("Final ground state");
    println("Up Density:");
    for(int j = 0; j < N; ++j)
        printfln("%d %.10f",1+j,upd(j));
    println();

    println("Dn Density:");
    for(int j = 0; j < N; ++j)
        printfln("%d %.10f",1+j,dnd(j));
    println();

    println("Total Density:");
    for(int j = 0; j < N; ++j)
        printfln("%d %.10f",1+j,(upd(j)+dnd(j)));
    println("Local Sz:");
    for(int j = 0; j < N; ++j)
        printfln("%d %.10f",1+j,szi(j));
    println();
}
void
KaicalculateAF(MPS const& psi,  SiteSet const& sites)
{
//  Calculate order of pararmeter Kai = Sumover i { S_i dot S_i+1 cross S_i+2 (start)}

   auto N = sites.N();
   auto psio = psi;
   auto Kaiq= 0.0+0.0_i;
	auto kai=std::vector<Complex>(N);
    FILE *kaif = fopen("kaiAF.txt", "w");
	FILE *kaisumf = fopen("kaiAFsum.txt", "w");
	


for(int i1 = 2; i1 <= (N-2); i1=i1+2)
{
        //Print(psi.A(i1-1));
        //Print(psi.A(i1));
        //Print(psi.A(i1+1));
        psio.position(i1);
       
        //ITensor Sx1 = 0.5*(sites.op("Sp",i1) + sites.op("Sm",i1));
        //ITensor Sy1 = (1/(2_i))*(sites.op("Sp",i1) - sites.op("Sm",i1));
        //ITensor Sz1 = sites.op("Sz",i1);

        ITensor Sx1 = sites.op("Sx",i1);
        ITensor Sy1 = sites.op("Sy",i1);
        ITensor Sz1 = sites.op("Sz",i1);
        ITensor Sx2 = sites.op("Sx",i1+1);
        ITensor Sy2 = sites.op("Sy",i1+1);
        ITensor Sz2 = sites.op("Sz",i1+1);
        ITensor Sx3 = sites.op("Sx",i1+2);
        ITensor Sy3 = sites.op("Sy",i1+2);
        ITensor Sz3 = sites.op("Sz",i1+2);
        auto B1 = Sx1*Sy2*Sz3 - Sx1*Sz2*Sy3 - Sy1*Sx2*Sz3 + Sy1*Sz2*Sx3 + Sz1*Sx2*Sy3 - Sz1*Sy2*Sx3;
        auto wf1 = psio.A(i1)*psio.A(i1+1)*psio.A(i1+2);
        // compute <wf1|B1|wf1>
        kai[i1-1]=(dag(prime(wf1,Site)) * B1 * wf1).cplx();
        Kaiq += kai[i1-1];
}
for(int i1 = 1; i1 <= (N-2); i1=i1+2)
{
        //Print(psi.A(i1-1));
        //Print(psi.A(i1));
        //Print(psi.A(i1+1));
        psio.position(i1);
       
        //ITensor Sx1 = 0.5*(sites.op("Sp",i1) + sites.op("Sm",i1));
        //ITensor Sy1 = (1/(2_i))*(sites.op("Sp",i1) - sites.op("Sm",i1));
        //ITensor Sz1 = sites.op("Sz",i1);

        ITensor Sx1 = sites.op("Sx",i1);
        ITensor Sy1 = sites.op("Sy",i1);
        ITensor Sz1 = sites.op("Sz",i1);
        ITensor Sx2 = sites.op("Sx",i1+1);
        ITensor Sy2 = sites.op("Sy",i1+1);
        ITensor Sz2 = sites.op("Sz",i1+1);
        ITensor Sx3 = sites.op("Sx",i1+2);
        ITensor Sy3 = sites.op("Sy",i1+2);
        ITensor Sz3 = sites.op("Sz",i1+2);
        auto B1 = -(Sx1*Sy2*Sz3 - Sx1*Sz2*Sy3 - Sy1*Sx2*Sz3 + Sy1*Sz2*Sx3 + Sz1*Sx2*Sy3 - Sz1*Sy2*Sx3);
        auto wf1 = psio.A(i1)*psio.A(i1+1)*psio.A(i1+2);
        // compute <wf1|B1|wf1>
        kai[i1-1]=(dag(prime(wf1,Site)) * B1 * wf1).cplx();
        Kaiq += kai[i1-1];
}
  printfln("\nsumKaiq = %.10f ",Kaiq);
  fprintf(kaisumf, "\nreal=%.10f \nimag=%.10f", Kaiq.real(),Kaiq.imag());
  println("\n ===============");
    for(int j = 1; j <= N; ++j)
    {
     printf("\n %.10f  ",kai[j-1].real());
	 fprintf(kaif, "%f\n", kai[j-1].real());
    }
  println("\n ===============");
  fclose(kaif);
  fclose(kaisumf);

}
MPO
makeMF1(Spinless const& sites)
    {
    auto N = sites.N();

    auto S2 = IQMPO(sites);

    auto links = std::vector<IQIndex>(N+1);
    for(auto n : range(N+1))
        {
        links.at(n) = IQIndex(nameint("L",n),
                              Index("0",3),QN("Nf=",0),
                              Index("+",1),QN("Nf=",0),
                              Index("-",1),QN("Nf=",0));
//							  Index("0",3),QN("Nf=",0,"Sz=",0),
//                              Index("+",1),QN("Nf=",0,"Sz=",-2),
//                              Index("-",1),QN("Nf=",0,"Sz=",+2));
        }

    for(auto n : range1(N))
        {
        auto row = dag(links.at(n-1));
        auto col = links.at(n);
        auto& W = S2.Aref(n);
        W = IQTensor(row,col,dag(sites(n)),prime(sites(n)));		
        W += sites.op("Id",n) * row(1) * col(1);
        W += sites.op("Id",n) * row(2) * col(2);
        W += sites.op("Id",n) * row(3) * col(3);
        W += sites.op("Id",n) * row(4) * col(4);
	W += sites.op("Id",n) * row(5) * col(5);
	if(n==1)
	{
        W = IQTensor(row,col,dag(sites(n)),prime(sites(n)));		
        W += sites.op("A",n) * row(1) * col(1);
        W += sites.op("Adag",n) * row(2) * col(2);
        W += sites.op("Id",n) * row(3) * col(3);
        W += sites.op("Id",n) * row(4) * col(4);
	W += sites.op("Id",n) * row(5) * col(5);
	}
  
/*
        W += sites.op("Id",n) * row(1) * col(1);
        W += sites.op("Id",n) * row(2) * col(2);

        W += sites.op("S2",n) * row(2) * col(1);

        W += 2*sites.op("Sz",n) * row(2) * col(3);
        W += sites.op("Id",n) * row(3) * col(3);
        W += sites.op("Sz",n) * row(3) * col(1);

        W += sites.op("S+",n) * row(2) * col(4);
        W += sites.op("Id",n) * row(4) * col(4);
        W += sites.op("S-",n) * row(4) * col(1);

        W += sites.op("S-",n) * row(2) * col(5);
        W += sites.op("Id",n) * row(5) * col(5);
        W += sites.op("S+",n) * row(5) * col(1);
*/
        W.scaleTo(1.);
        }

 //   S2.Aref(1) *= setElt(links.at(0)(2)).;
    auto rindex=links.at(0);
    auto lindex=dag(links.at(N));
    auto rvector=IQTensor(rindex);
    auto lvector=IQTensor(lindex);
    rvector.set(rindex(1),1.0);
    rvector.set(rindex(2),1.0);
    lvector.set(lindex(1),1.0);
    lvector.set(lindex(2),1.0);


    S2.Aref(1) *= rvector;
    S2.Aref(N) *= lvector;
   auto S2MPO = MPO(sites);
    for(int j = 1; j <= N; ++j)
	{
	auto A = S2.A(j);
	S2MPO.setA(j,toITensor(A));
	}
    return S2MPO;
  //  return S2;
    }
MPO
makeMF2(Spinless const& sites)
    {
    auto N = sites.N();

    auto S2 = IQMPO(sites);

    auto links = std::vector<IQIndex>(N+1);
    for(auto n : range(N+1))
        {
        links.at(n) = IQIndex(nameint("L",n),
                              Index("0",3),QN("Nf=",0),
                              Index("+",1),QN("Nf=",0),
                              Index("-",1),QN("Nf=",0));
//							  Index("0",3),QN("Nf=",0,"Sz=",0),
//                              Index("+",1),QN("Nf=",0,"Sz=",-2),
//                              Index("-",1),QN("Nf=",0,"Sz=",+2));
        }

    for(auto n : range1(N))
        {
        auto row = dag(links.at(n-1));
        auto col = links.at(n);
        auto& W = S2.Aref(n);
        W = IQTensor(row,col,dag(sites(n)),prime(sites(n)));		
        W += sites.op("Id",n) * row(1) * col(1);
        W += sites.op("Id",n) * row(2) * col(2);
        W += sites.op("Id",n) * row(3) * col(3);
        W += sites.op("Id",n) * row(4) * col(4);
	W += sites.op("Id",n) * row(5) * col(5);
	if(n==1)
	{
        W = IQTensor(row,col,dag(sites(n)),prime(sites(n)));		
        W += (1.0/1.0_i)*sites.op("A",n) * row(1) * col(1);
        W += (-1.0/1.0_i)*sites.op("Adag",n) * row(2) * col(2);
        W += sites.op("Id",n) * row(3) * col(3);
        W += sites.op("Id",n) * row(4) * col(4);
	W += sites.op("Id",n) * row(5) * col(5);
	}
  
/*
        W += sites.op("Id",n) * row(1) * col(1);
        W += sites.op("Id",n) * row(2) * col(2);

        W += sites.op("S2",n) * row(2) * col(1);

        W += 2*sites.op("Sz",n) * row(2) * col(3);
        W += sites.op("Id",n) * row(3) * col(3);
        W += sites.op("Sz",n) * row(3) * col(1);

        W += sites.op("S+",n) * row(2) * col(4);
        W += sites.op("Id",n) * row(4) * col(4);
        W += sites.op("S-",n) * row(4) * col(1);

        W += sites.op("S-",n) * row(2) * col(5);
        W += sites.op("Id",n) * row(5) * col(5);
        W += sites.op("S+",n) * row(5) * col(1);
*/
        W.scaleTo(1.);
        }

 //   S2.Aref(1) *= setElt(links.at(0)(2)).;
    auto rindex=links.at(0);
    auto lindex=dag(links.at(N));
    auto rvector=IQTensor(rindex);
    auto lvector=IQTensor(lindex);
    rvector.set(rindex(1),1.0);
    rvector.set(rindex(2),1.0);
    lvector.set(lindex(1),1.0);
    lvector.set(lindex(2),1.0);


    S2.Aref(1) *= rvector;
    S2.Aref(N) *= lvector;
   auto S2MPO = MPO(sites);
    for(int j = 1; j <= N; ++j)
	{
	auto A = S2.A(j);
	S2MPO.setA(j,toITensor(A));
	}
    return S2MPO;
  //  return S2;
    }

void
printmatrixstate(MPS const& psiq, int N)
{
   printfln("matrixstate\n");
   ITensor phiq=psiq.A(1);
     for(int j = 2; j <= N; ++j)
        {
         phiq = phiq*(psiq.A(j));
        }

    PrintData(phiq);
}
void
printmatrixmpo(MPO const& psiq, int N)
{
   printfln("matrixHamltonian\n");
   ITensor phiq=psiq.A(1);
     for(int j = 2; j <= N; ++j)
        {
         phiq = phiq*(psiq.A(j));
        }

    PrintData(phiq);
}
void
measureobs(MPS psi,SiteSet const& sites)
{ 
    auto N=sites.N();
    auto Nsite=std::vector<Complex>(N);
    auto mf1d=std::vector<Complex>(N);
    auto mf2d=std::vector<Complex>(N);
     auto exA = 0.0+0.0_i;   
     auto exAdag = 0.0+0.0_i;   
     auto exN = 0.0+0.0_i;   
    for(int j = 1; j <= N; ++j)
        {
        psi.position(j);
        exA = (dag(prime(psi.A(j),Site))*sites.op("C",j)*psi.A(j)).cplx();
	exAdag = (dag(prime(psi.A(j),Site))*sites.op("Cdag",j)*psi.A(j)).cplx();
        exN = (dag(prime(psi.A(j),Site))*sites.op("N",j)*psi.A(j)).cplx();
        mf1d[j-1] = exA + exAdag;
        mf2d[j-1] = 1.0_i*exAdag-1.0_i*exA;
        Nsite[j-1] = exN;
        }


    println("Ni:");
    for(int j = 0; j < N; ++j)
        printfln("%d %.10f",1+j,Nsite[j]);
    println();
}
MPO
makeS2(Hubbard const& sites)
    {
    auto N = sites.N();

    auto S2 = IQMPO(sites);

    auto links = std::vector<IQIndex>(N+1);
    for(auto n : range(N+1))
        {
        links.at(n) = IQIndex(nameint("L",n),
                              Index("0",3),QN("Nf=",0,"Sz=",0),
                              Index("+",1),QN("Nf=",0,"Sz=",-2),
                              Index("-",1),QN("Nf=",0,"Sz=",+2));
        }

    for(auto n : range1(N))
        {
        auto row = dag(links.at(n-1));
        auto col = links.at(n);
        auto& W = S2.Aref(n);
        W = IQTensor(row,col,dag(sites(n)),prime(sites(n)));

        W += sites.op("Id",n) * row(1) * col(1);
        W += sites.op("Id",n) * row(2) * col(2);

        W += sites.op("S2",n) * row(2) * col(1);

        W += 2*sites.op("Sz",n) * row(2) * col(3);
        W += sites.op("Id",n) * row(3) * col(3);
        W += sites.op("Sz",n) * row(3) * col(1);

        W += sites.op("S+",n) * row(2) * col(4);
        W += sites.op("Id",n) * row(4) * col(4);
        W += sites.op("S-",n) * row(4) * col(1);

        W += sites.op("S-",n) * row(2) * col(5);
        W += sites.op("Id",n) * row(5) * col(5);
        W += sites.op("S+",n) * row(5) * col(1);

        W.scaleTo(1.);
        }

    S2.Aref(1) *= setElt(links.at(0)(2));
    S2.Aref(N) *= setElt(dag(links.at(N))(1));
    auto S2MPO = MPO(sites);
    for(int j = 1; j <= N; ++j)
        {
	  auto A = S2.A(j);
	  S2MPO.setA(j,toITensor(A));
	      }
    return S2MPO;
}
void
caledgeSzSz(Hubbard const& sites,MPS psi)
{
  //constructed from a SiteSet "sites"

//Replace "Op1" and "Op2" with the actual names
//of the operators you want to measure
auto N=sites.N();
printfln("\n the Sz-Sz correlation from floor(N/4)... floor(3N/4)  are: \n  ");
for(int j=int(N/4)+1; j<=int(N*3/4) ;++j)
{
  int i=N/4;
//  auto j=sites.N();
  auto op_i = sites.op("Sz",i);
  auto op_j = sites.op("Sz",j);

  
//below we will assume j > i

//'gauge' the MPS to site i
//any 'position' between i and j, inclusive, would work here
  psi.position(i); 

//psi.Anc(1) *= psi.A(0); //Uncomment if doing iDMRG calculation

//index linking i to i+1:
  auto ir = commonIndex(psi.A(i),psi.A(i+1),Link);

  auto C = psi.A(i)*op_i*dag(prime(psi.A(i),Site,ir));
  for(int k = i+1; k < j; ++k)
      {
      C *= psi.A(k);
      C *= dag(prime(psi.A(k),Link));
      }
    C *= psi.A(j);
    C *= op_j;
//index linking j to j-1:
    auto jl = commonIndex(psi.A(j),psi.A(j-1),Link);
    C *= dag(prime(psi.A(j),jl,Site));

    auto result = C.real(); //or C.cplx() if expecting complex
    printfln("%.10f",result);
}
}
void
calstringorder(Hubbard const& sites,MPS psi)
{
  //constructed from a SiteSet "sites"

//Replace "Op1" and "Op2" with the actual names
//of the operators you want to measure
  auto N=sites.N();
  char filename[20];
  FILE* fptr2;
  sprintf(filename,"string%d_%d.txt",detail::quickran(),N);
  fptr2 = fopen(filename,"w");


  printfln("\n the string order from 2 to N  are: \n ");
 // Vector vecstring(N-1);
  
for(int j=2; j<=N ;++j)
{
  int i=1;
//  auto j=sites.N();
  auto op_i = sites.op("Sz",i);
  auto op_j = sites.op("Sz",j);

  
//below we will assume j > i

//'gauge' the MPS to site i
//any 'position' between i and j, inclusive, would work here
  psi.position(i); 

//psi.Anc(1) *= psi.A(0); //Uncomment if doing iDMRG calculation

//index linking i to i+1:
  auto ir = commonIndex(psi.A(i),psi.A(i+1),Link);

  auto C = psi.A(i)*op_i*dag(prime(psi.A(i),Site,ir));
  for(int k = i+1; k < j; ++k)
      {
      C *= psi.A(k);
//      C *= dag(prime(psi.A(k),Link));
      C *= dag(prime(psi.A(k),Link,Site));
	    C *= (-1.0)*sites.op("F",k);
      }
    C *= psi.A(j);
    C *= op_j;
//index linking j to j-1:
    auto jl = commonIndex(psi.A(j),psi.A(j-1),Link);
    C *= dag(prime(psi.A(j),jl,Site));

    auto result = C.real(); //or C.cplx() if expecting complex
	fprintf(fptr2, "%.10f\n", result);

}
  printfln("\n check point before return vecstring\n ");
    fclose(fptr2);
}
double
enspec(Hubbard const& sites,MPS psi,int b)
{
    //auto b=sites.N()/2;
    psi.position(b); 

//Here assuming an MPS of ITensors, but same code works
//for IQMPS by replacing ITensor -> IQTensor

//Compute two-site wavefunction for sites (b,b+1)
    ITensor wf = psi.A(b)*psi.A(b+1);

//SVD this wavefunction to get the spectrum
//of density-matrix eigenvalues
    auto UL = psi.A(b);
    ITensor S,VR;
    auto spectrum = svd(wf,UL,S,VR);
//Apply von Neumann formula
//spectrum.eigs() is a Vector containing
//the density matrix eigenvalues
//(squares of the singular values)
    Real SvN = 0.;
//    printfln("Entanglement spectrum");
    for(auto p : spectrum.eigs())
        {
//        printfln("%.10f",p);
        if(p > 1E-8) 
			{
//				printfln("%.10f",p);
				SvN += -p*log(p);
			}
        }
 //   printfln("Across bond b=%d, SvN = %.10f",b,SvN);
	return SvN;
}
int main(int argc, char* argv[])
    {
    //Parse the input file
    if(argc != 2) { printfln("Usage: %s inputfile_exthubbard",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");

    auto N = input.getInt("N");
    auto NEex = input.getInt("NEex");
    auto Npart = input.getInt("Npart",N); //number of particles, default is N (half filling)
    auto wei1 = input.getReal("wei1",20.0);
    auto nsweeps = input.getInt("nsweeps");
    auto t1 = input.getReal("t1",1);
    auto t2 = input.getReal("t2",1);
    auto t3 = input.getReal("t3",1);
    auto t4 = input.getReal("t4",1);
    auto bz = input.getReal("bz",0);
    auto tso = input.getReal("tso",0);
    auto pinNup = input.getReal("pinNup",0);
    auto pinNdn = input.getReal("pinNdn",0);
    auto pinSz = input.getReal("pinSz",0);
    auto pinSx = input.getReal("pinSx",0);
    auto pinSy = input.getReal("pinSy",0);
    auto U = input.getReal("U",0);
    auto Uo = input.getReal("Uo",0);
    auto Vex = input.getReal("Vex",0);
    auto mu = input.getReal("mu",0);
    auto Jz = input.getReal("Jz",0);
    auto Jxx = input.getReal("Jxx",0);
    auto hx = input.getReal("hx",0);
    auto hz = input.getReal("hz",0);
	auto hzr = input.getReal("hzr",0);
	auto hxr = input.getReal("hxr",0);
	auto nr = input.getReal("nr",0);		
    auto quiet = input.getYesNo("quiet",false);
    auto writestates=input.getYesNo("writestates",false);
    auto table = InputGroup(input,"sweeps");
    auto sweeps = Sweeps(nsweeps,table);

    println(sweeps);

    //
    // Initialize the site degrees of freedom.
    //
  //  auto sites = Hubbard(N);
//    auto sites =Hubbard(N,{"ConserveSz",true});
    auto sites = Hubbard(N,{"ConserveQNs",false});
    //
    // Create the Hamiltonian using AutoMPO
    //
    auto ampo = AutoMPO(sites);
 //   for(int i = 1; i <= N; ++i) 
 //       {
  //      ampo += U,"Nupdn",i;
   //     }
	Vector himp(N),nimp(N);
	 for(int i = 0; i < N; i=i+1)
	 {
		 himp(i)=detail::quickran()*2.0-1.0;// in range -1 to 1
		 nimp(i)=detail::quickran()*2.0-1.0;
//		 printfln("%.10f",himp(i));
//		 printfln("%.10f",nimp(i));
	 }
    for(int b = 1; b <= (N-1); b=b+1)// for chain#1,
        {
       ampo += -t1,"Cdagup",b,"Cup",b+1;
        ampo += -t1,"Cdagup",b+1,"Cup",b; //cc
        ampo += -t1,"Cdagdn",b,"Cdn",b+1;
        ampo += -t1,"Cdagdn",b+1,"Cdn",b;//cc
 //       ampo += V1,"Ntot",b,"Ntot",b+1;
        }


            println("\nchecking point after t1 ");



 

   

////////////////////////////////////
    for(int i = 1; i <=(N); i=i+1)// imp spin sz
       {
//		printfln("%.10f", hzr*himp(i));
        ampo += hzr*himp(i-1),"Sz",i;
        }
//////////////////////////////////////////imp spin Sx		
    for(int i = 1; i <=(N); i=i+1)// 
       {
        ampo += (hxr*himp(i-1))/2.0,"S+",i;
        }
    for(int i = 1; i <=(N); i=i+1)// 
       {
        ampo += (hxr*himp(i-1))/2.0,"S-",i;
        }
//////////////////////////////////////////////
////////////////////////////////////
    for(int i = 1; i <=(N); i=i+1)// imp density onsite
       {
        ampo += nr*nimp(i-1),"Ntot",i;
        }
///////////////////////////////////////////////
    for(int i = 1; i <=(N); i=i+1)// on site U
       {
        ampo += U,"Nup",i,"Ndn",i;
        }
    println("\nchecking point after U ");
    for(int i = 1; i <=(N); i=i+1)// -mu*N chemical potential
       {
        ampo += -mu,"Ntot",i;
        }
    println("\nchecking point after mu ");

/////////////////////////////////////
    for(int i = 1; i <=(N-1); i=i+1)// on site U
       {
		ampo += (0.5)*Jxx,"S+",i,"S-",i+1;
        ampo += (0.5)*Jxx,"S-",i,"S+",i+1;
        ampo += Jz,"Sz",i,"Sz",i+1;
        }
    for(int i = 1; i <=(N); i=i+1)// on site U
       {
        ampo += hx/2.0,"S+",i;
        }
    for(int i = 1; i <=(N); i=i+1)// on site U
       {
        ampo += hx/2.0,"S-",i;
        }
    for(int i = 1; i <=(N); i=i+1)// on site U
       {
        ampo += hz,"Sz",i;
        }
////////////////////////////////

		ampo += pinNup,"Nup",1;
		ampo += pinNdn,"Ndn",1;
		ampo += pinSz,"Sz",1;

		ampo +=0.5*pinSx,"S+",1;
		ampo += 0.5*pinSx,"S-",1;
		ampo +=0.5_i* pinSy,"S+",1;
		ampo += -0.5_i*pinSy,"S-",1;


//        ampo += pin,"N",2;
//        ampo += pin,"N",3;
//        ampo += pin,"N",4;
            println("\nchecking point before IQMPO ");
    auto H = MPO(ampo);

            println("\nchecking point after IQMPO ");

    //
    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    auto state = InitState(sites);
            println("\nchecking point after state declar ");
    int p = Npart;
    for(int i = N; i >= 1; --i) 
        {
        if(p > i)
            {
            println("Doubly occupying site ",i);
            state.set(i,"UpDn");
            p -= 2;
            }
        else
        if(p > 0)
            {
            println("Singly occupying site ",i);
            state.set(i,(i%2==1 ? "Up" : "Dn"));
            p -= 1;
            }
        else
            {
            state.set(i,"Emp");
            }
        }


    // Begin the DMRG calculation
//	auto psi = MPS(state);
	auto psid = std::vector<MPS>(50);
    for(int d=0;d<=49;d++) psid.at(d) = MPS(state);
	
    Vector envec(NEex+1);
		/*
	Matrix enprofile(NEex+1,N-1);

	Matrix stringorder(NEex+1,N-1);
	Vector stringvec(N-1);
	*/
    envec(0) = dmrg(psid.at(0),H,sweeps,{"Quiet",quiet});
	/*	
	for(int i=1; i<(N);i++)
	{
			enprofile(0,i-1)=enspec(sites,psid.at(0),i);
	}

	stringvec=calstringorder(sites,psid.at(0));
	for(int i=0; i<(N);i++)
	{
			stringorder(0,i)=stringvec(i);
			
	}
	printfln("\n check point after ground state cal\n ");
	
	*/

   for(int d=1;d<=NEex;d++)
   {
	
    auto wfsqd = std::vector<MPS>(d);
	printfln("\n check point after declared wfsqd\n ");
    for(int i=d;i>0;i--) wfsqd.at(d-i)= psid.at(d-i);    
	printfln("\n check point after filling wfsqd with psid\n ");
	printfln("\n check point print size of wfsqd\n ",wfsqd.size());
    envec(d) = dmrg(psid.at(d),H,wfsqd,sweeps,{"Quiet",true,"Weight",wei1});
    println("entanlement");
	/*
	for(int i=1; i<N;i++)
	{
			enprofile(d,i-1)=enspec(sites,psid.at(d),i);
	}
	
	stringvec=calstringorder(sites,psid.at(d));
	for(int i=0; i<(N);i++)
	{
			stringorder(d,i)=stringvec(i);
	}
 
	*/
	
    if(writestates==true)
    {
    writeToFile(format("sites_%d",N),sites); //file name will be sites_100
      for(int k=0;k<d;++k)
      {
      writeToFile(format("psi_%d",k),psid.at(k));
      }
    } 	
   
    for(int i=(d-1);i>=0;i--)
    {    
    auto otho20 = overlapC(psid.at(d),psid.at(i));//check the orthognality
    printf("\n |<psiq%d|psiq%d>|=%.10f",i,d,abs(otho20));
    }
   }

  
	std::vector<std::pair<double, int> > vp; 
	for (int i = 0; i < NEex+1; ++i) { 
        vp.push_back(std::make_pair(envec(i), i)); 
    }
	
//    std::sort(envec.begin(),envec.end());
    sort(vp.begin(), vp.end()); 
	auto Eoindex=vp[0].second;
	printf("\n Eo = %.10f index=%d",vp[0].first,vp[0].second);
	println();
    auto gap=vp[1].first-vp[0].first;
    printfln("gap = %.10f",gap);  
	
	
	char filename[20];
	FILE* fptr;

	sprintf(filename,"Svnprofile%.10f_%d.txt",detail::quickran(),N);
	fptr = fopen(filename,"w");


	
	for(int i=1; i<(N);i++)
	{
			fprintf(fptr, "%.10f\n", enspec(sites,psid.at(Eoindex),i));
	}

	calstringorder(sites,psid.at(Eoindex));
	
	fclose(fptr);

	
	
    for(int j = 0; j < (NEex+1); ++j) printfln("%.10f",envec(j));






    return 0;
    }
