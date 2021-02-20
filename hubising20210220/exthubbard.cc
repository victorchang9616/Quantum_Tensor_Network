#include "itensor/all.h"
#include "itensor/util/stdx.h"
using namespace itensor;
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
  printfln("\n the string order from floor(N/4)... floor(3N/4)  are: \n ");
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
    printfln("%.10f",result);
}
}
void
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
				printfln("%.10f",p);
				SvN += -p*log(p);
			}
        }
    printfln("Across bond b=%d, SvN = %.10f",b,SvN);
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

            println("\nchecking point after stateassign ");
    auto psi = MPS(state);

            println("\nchecking point after IQMPS ");
    Vector upd(N),dnd(N),szi(N);
/*
    for(int j = 1; j <= N; ++j)
        {
        psi.position(j);
        upd(j-1) = (dag(prime(psi.A(j),Site))*sites.op("Nup",j)*psi.A(j)).real();
        dnd(j-1) = (dag(prime(psi.A(j),Site))*sites.op("Ndn",j)*psi.A(j)).real();
        szi(j-1) = (dag(prime(psi.A(j),Site))*sites.op("Sz",j)*psi.A(j)).real();
        }

    println("Initial state");
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
*/
//    Print(totalQN(psi));

    //
    // Begin the DMRG calculation
    Vector envec(NEex+1);
    envec(0) = dmrg(psi,H,sweeps,{"Quiet",quiet,"WriteM",4000});
	    printfln("Entanglement spectrum");
	for(int i=1; i<(N);i++)
	{
			enspec(sites,psi,i);
	}

 //   Print(totalQN(psi));
    //
    // Measure spin densities
    //
  //  Vector upd(N),dnd(N),szi(N);
	caledgeSzSz(sites,psi);
//    printfln("Entanglement spectrum");
	calstringorder(sites,psi);

    auto S2MPO = makeS2(sites);
    auto S20 = overlap(psi,S2MPO,psi);
    auto S0=(-1.0+sqrt(1+4*S20))/2.0;
    printfln("\nS2 = %.10f",S20);
	  printfln("\nS = %.10f",S0);
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
//and some particular bond "b" (1 <= b < psi.N())
//across which we want to compute the von Neumann entanglement

//	enspec(sites,psi);
    //
    // Print the final energy reported by DMRG
    //
//    printfln("\nGround State Energy = %.10f",energy);
	auto psi1 = MPS(state);
	auto wfsq = std::vector<MPS>(1);
    wfsq.at(0) = psi;
   envec(1)  = dmrg(psi1,H,wfsq,sweeps,{"Quiet",true,"Weight",wei1,"WriteM",4000});
   // measurement
   println("entanlement");
	for(int i=1; i<N;i++)
	{
			enspec(sites,psi1,i);
	}
   caledgeSzSz(sites,psi1);
   calstringorder(sites,psi1);
   S20 = overlap(psi1,S2MPO,psi1);
   S0=(-1.0+sqrt(1+4*S20))/2.0;
   printfln("\nS2 = %.10f",S20);
   printfln("\nS = %.10f",S0);
   for(int j = 1; j <= N; ++j)
        {
        psi1.position(j);
        upd(j-1) = (dag(prime(psi1.A(j),Site))*sites.op("Nup",j)*psi1.A(j)).real();
        dnd(j-1) = (dag(prime(psi1.A(j),Site))*sites.op("Ndn",j)*psi1.A(j)).real();
        szi(j-1) = (dag(prime(psi1.A(j),Site))*sites.op("Sz",j)*psi1.A(j)).real();
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

   // measurement
 //   writeToFile(format("psi1_%d",N),psi1);
//    auto psi1 = iqmpstomps(psiq1,sites);
    auto wfsq2 = std::vector<MPS>(2);
    wfsq2.at(0) = psi;
    wfsq2.at(1) = psi1;
 //   auto S21 = overlap(psi0,S2MPO,psi0);
  //  auto S1=(-1.0+sqrt(1+4*S21))/2.0;
  //  printfln("\nen1 = %.10f, S2 = %.10f",en1,S21);
//	printfln("\nS = %.10f",S1);
 //   Kaicalculate(psi1,sites);
//    trimerbond(K11,K12,K21,K22,psi0,sites);
   // Sicalculate(psi1,sites);
    auto otho = overlapC(psi,psi1);
    println("\n |<psiq0|psiq1>|",abs(otho));
       auto psid = std::vector<MPS>(50);
    psid.at(0) = psi;
    psid.at(1) = psi1; 
    for(int d=2;d<=49;d++) psid.at(d) = MPS(state);

   for(int d=2;d<=NEex;d++)
   {

    auto wfsqd = std::vector<MPS>(d);
    for(int i=d;i>0;i--) wfsqd.at(d-i)= psid.at(d-i);    
    envec(d) = dmrg(psid.at(d),H,wfsqd,sweeps,{"Quiet",true,"Weight",wei1,"WriteM",4000});
    println("entanlement");
	for(int i=1; i<N;i++)
	{
			enspec(sites,psid.at(d),i);
	}
	calstringorder(sites,psid.at(d));
	caledgeSzSz(sites,psid.at(d));
    S20 = overlap(psid.at(d),S2MPO,psid.at(d));
    S0=(-1.0+sqrt(1+4*S20))/2.0;
    printfln("\nS2 = %.10f",S20);
    printfln("\nS = %.10f",S0);   
 // measurement
      for(int j = 1; j <= N; ++j)
        {
        psid.at(d).position(j);
        upd(j-1) = (dag(prime(psid.at(d).A(j),Site))*sites.op("Nup",j)*psid.at(d).A(j)).real();
        dnd(j-1) = (dag(prime(psid.at(d).A(j),Site))*sites.op("Ndn",j)*psid.at(d).A(j)).real();
        szi(j-1) = (dag(prime(psid.at(d).A(j),Site))*sites.op("Sz",j)*psid.at(d).A(j)).real();
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
	
    if(writestates==true)
    {
    writeToFile(format("sites_%d",N),sites); //file name will be sites_100
      for(int k=0;k<d;++k)
      {
      writeToFile(format("psi_%d",k),psid.at(k));
      }
    } 	
  //  writeToFile(format("psi%d_%d",d,N),psid.at(d));
//    auto psi1 = iqmpstomps(psiq1,sites);


  //  auto S22 = overlap(psid.at(d),S2MPO,psid.at(d));
   // auto S2=(-1.0+sqrt(1+4*S22))/2.0;
   // printfln("\nen%d = %.10f, S2 = %.10f",d,en2,S22);
   // printfln("\nS = %.10f",S2);
   // Kaicalculate(psid.at(d),sites);
//    trimerbond(K11,K12,K21,K22,psi0,sites);
   // Sicalculate(psid.at(d),sites);
   
    for(int i=(d-1);i>=0;i--)
    {    
    auto otho20 = overlapC(psid.at(d),psid.at(i));//check the orthognality
    printf("\n |<psiq%d|psiq%d>|=%.10f",i,d,abs(otho20));
    }
   }
    std::sort(envec.begin(),envec.end());
	println();
    auto gap=envec(1)-envec(0);
    printfln("gap = %.10f",gap);
    for(int j = 0; j < (NEex+1); ++j)
	printfln("%.10f",envec(j));
    return 0;
    }
