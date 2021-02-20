#include "itensor/all.h"
#include <list>  
#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651
using namespace itensor;
//using std::vector;
using namespace std; 

void 
Swap(Real array[], int x, int y)
{
    int temp = array[x];
    array[x] = array[y];
    array[y] = temp;
}

void
BubbleSort(Real array[], int size)
{
  for(int i=0; i < size; i++)
  {
    for(int j=1; j < size - i; j++)
    {
       if(array[j] < array[j-1])
       {
          Swap(array, j, j-1);
       }
    }
  }
}
IQMPS
Singlet(SiteSet const& sites)
{

 // singlets for 2N sites 2N -> 4 here
    auto psis = IQMPS(sites);
    for(int n = 1; n <= 4; n += 2)
    {
    auto s1 = sites(n);
    auto s2 = sites(n+1);
    auto wf = IQTensor(s1,s2);
    wf.set(s1(1),s2(2), ISqrt2);
    wf.set(s1(2),s2(1), -ISqrt2);
    IQTensor D;
    psis.Aref(n) = IQTensor(s1);
    psis.Aref(n+1) = IQTensor(s2);
    svd(wf,psis.Aref(n),D,psis.Aref(n+1));
    psis.Aref(n) *= D;
    }
    return psis;
}

ITensor
IQMPOtoITensor(SiteSet const& sites,IQMPO const& Hq)
{
	auto N = sites.N();
	IQTensor Hqa = Hq.A(1);
	for(int j = 2; j <= N; ++j)
	{
         Hqa = Hqa*(Hq.A(j));
	}
	auto Ha = toITensor(Hqa);
	return Ha;
}
ITensor
MPOtoITensor(SiteSet const& sites,MPO const& Hq)
{
	auto N = sites.N();
	ITensor Hqa = Hq.A(1);
	for(int j = 2; j <= N; ++j)
	{
         Hqa = Hqa*(Hq.A(j));
	}
	auto Ha = Hqa;
	return Ha;
}
IQMPO
makeSz(SiteSet const& sites)
{
    auto N = sites.N();
    auto ampo = AutoMPO(sites);// for H1D
   // auto ampo1 = AutoMPO(sites);//for Hp

      for(int j = 1; j <= (N); ++j)
        {
         ampo += 1.0,"Sz",j;
        }
     auto Sz = IQMPO(ampo);
     return Sz;     
}
IQMPO
makeSz2(SiteSet const& sites)
{
    auto N = sites.N();
    auto ampo = AutoMPO(sites);// for H1D
   // auto ampo1 = AutoMPO(sites);//for Hp

      for(int j = 1; j <= (N); ++j)
        {
         ampo += 1.0,"Sz*Sz",j;
        }
     auto Sz2 = IQMPO(ampo);
     return Sz2;     
}
IQMPO
makeSx2(SiteSet const& sites)
{
    auto N = sites.N();
    auto ampo = AutoMPO(sites);// for H1D
   // auto ampo1 = AutoMPO(sites);//for Hp

      for(int j = 1; j <= (N); ++j)
        {
         ampo += 1.0,"Sp*Sp",j;
         ampo += 1.0,"Sm*Sm",j;
         ampo += 1.0,"Sp*Sm",j;
         ampo += 1.0,"Sm*Sp",j;
        }
     auto Sx2 = IQMPO(ampo);
     return Sx2;     
}
IQMPO
makeSy2(SiteSet const& sites)
{
    auto N = sites.N();
    auto ampo = AutoMPO(sites);// for H1D
   // auto ampo1 = AutoMPO(sites);//for Hp

      for(int j = 1; j <= (N); ++j)
        {
         ampo += 1.0,"Sy*Sy",j;
        }
     auto Sy2 = IQMPO(ampo);
     return Sy2;     
}
IQMPO
makekai(SiteSet const& sites)
{
    auto N = sites.N();
    auto ampo = AutoMPO(sites);// for H1D
   // auto ampo1 = AutoMPO(sites);//for Hp

      for(int j = 1; j <= (N-2); ++j)
        {
         ampo += (0.5_i),"Sz",j,"Sp",j+1,"Sm",j+2;
         ampo += -(0.5_i),"Sz",j,"Sm",j+1,"Sp",j+2;
         ampo += -(0.5_i),"Sp",j,"Sz",j+1,"Sm",j+2;
         ampo += (0.5_i),"Sm",j,"Sz",j+1,"Sp",j+2;
         ampo += (0.5_i),"Sp",j,"Sm",j+1,"Sz",j+2;
         ampo += -(0.5_i),"Sm",j,"Sp",j+1,"Sz",j+2;

        }
     auto kai = IQMPO(ampo);
     return kai;     
}
IQMPO
makeS2(SiteSet const& sites)
    {
    auto N = sites.N();

    auto S2 = IQMPO(sites);

    auto links = std::vector<IQIndex>(N+1);
    for(auto n : range(N+1))
        {
 
        links.at(n) = IQIndex(nameint("L",n),
                              Index("0",3),QN("Sz=",0),
                              Index("+",1),QN("Sz=",-2),
                              Index("-",1),QN("Sz=",+2));
 /*       links.at(n) = IQIndex(nameint("L",n),
                              Index("0",3),QN("Nf=",0,"Sz=",0),
                              Index("+",1),QN("Nf=",0,"Sz=",-2),
                              Index("-",1),QN("Nf=",0,"Sz=",+2));
*/
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

    return S2;
    }

MPO
makeS2mpo(SiteSet const& sites)
    {
    auto N = sites.N();

    auto S2 = IQMPO(sites);

    auto links = std::vector<IQIndex>(N+1);
    for(auto n : range(N+1))
        {

        links.at(n) = IQIndex(nameint("L",n),
                              Index("0",3),QN("Sz=",0),
                              Index("+",1),QN("Sz=",-2),
                              Index("-",1),QN("Sz=",+2));
/*
        links.at(n) = IQIndex(nameint("L",n),
                              Index("0",3),QN("Nf=",0,"Sz=",0),
                              Index("+",1),QN("Nf=",0,"Sz=",-2),
                              Index("-",1),QN("Nf=",0,"Sz=",+2));
*/
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
Sicalculate(IQMPS const& psi, SiteSet const& sites)
{
// Compute magnetization for ITensor
    auto N = sites.N();
    auto psio = psi;
    auto Sxj = 0.0+0.0_i;
    auto Syj = 0.0+0.0_i; 
    auto Szj = 0.0+0.0_i;
    auto S2j = 0.0+0.0_i;
    for(int j = 1; j <= N; ++j)
    {
     ITensor Sx = sites.op("Sx",j);
     ITensor Sy = sites.op("Sy",j);
     ITensor Sz = sites.op("Sz",j);
     ITensor S2 = sites.op("S2",j);
    //Make site j the MPS "orthogonality center"
    psio.position(j);
    //Measure magnetization
     auto Sxex = (toITensor(psio.A(j))* Sx * dag(prime(toITensor(psio.A(j)),Site))).cplx();
     auto Syex = (toITensor(psio.A(j))* Sy * dag(prime(toITensor(psio.A(j)),Site))).cplx();
     auto Szex = (toITensor(psio.A(j))* Sz * dag(prime(toITensor(psio.A(j)),Site))).cplx();
     auto S2ex = (toITensor(psio.A(j))* S2 * dag(prime(toITensor(psio.A(j)),Site))).cplx();

    Sxj = Sxj + Sxex;
    Syj = Syj + Syex;
    Szj = Szj + Szex;
    S2j = S2j + S2ex;
   // Syj = Syj + (psi.A(j)* Sy * dag(prime(psi.A(j),Site))).cplx();
    //Szj = Szj + (psi.A(j)* sites.op("Sz",j)* dag(prime(psi.A(j),Site))).cplx();
   //S2j = S2j + (psi.A(j)* S2 * dag(prime(psi.A(j),Site))).cplx();
    //println("Sz_",j," = ",Szj);
    printf("\nSx%d = %.5f  ",j,Sxex);
    printf("\nSy%d = %.5f",j,Syex);
    printf("\nSz%d = %.5f",j,Szex);
    printf("\nS2%d = %.5f ",j,S2ex);
    }
   println("\n===============");
    println("\nSxtotal =  ",Sxj);
    println("\nSytotal =  ",Syj);
    println("\nSztotal =  ",Szj);
    println("\nS2total =  ",S2j);
    println("\n ===============");
}
void
Kaicalculate(IQMPS const& psi,  SiteSet const& sites)
{
//  Calculate order of pararmeter Kai = Sumover i { S_i dot S_i+1 cross S_i+2 (start)}
   auto N = sites.N();
   auto psio = psi;
   auto Kaiq= 0.0+0.0_i;

for(int i1 = 1; i1 <= (N-2); ++i1)
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
        auto wf1 = toITensor(psio.A(i1))*toITensor(psio.A(i1+1))*toITensor(psio.A(i1+2));
        // compute <wf1|B1|wf1>
        Kaiq +=  (dag(prime(wf1,Site)) * B1 * wf1).cplx();
}
  println("\nKaiq = ",Kaiq);
}
/*
void
S2calculate(IQMPS const& psi, SiteSet const& sites)
{
  auto N = sites.N();
  ITensor Szsum = sites.op("Sz",1);
  for(int i1 = 2; i1 <= N; ++i1)
  {
	Szsum = Szsum + sites.op("Sz",i1)
  }
  PrintData(Szsum);
}
*/
void
printmatrixstateq(IQMPS const& psiq, int N)
{
   printfln("IQmatrixstate\n");
   IQTensor phiq=psiq.A(1);
     for(int j = 2; j <= N; ++j)
        {
         phiq = phiq*(psiq.A(j));
        }

    PrintData(phiq);
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
timereversal(MPS const& psi,SiteSet const& sites)
{
    // check TRS of psi////////////////////////////////////
    auto N = sites.N();
    MPS psit = psi;
    auto vinds = std::vector<Index>(3);
for(int k=1; k <= N; ++k)
{
    auto T= psit.A(k); // Site n
    int iind=1;
    for(auto& I : T.inds())
{   

   if(I.type()==Link) 
    {
     vinds[iind] = I;
     iind++;
    }
    if(I.type()==Site) vinds[0]= I;
}
 
    auto St=vinds[0];
    auto b1=vinds[1];
    auto b2=vinds[2];
if(T.r() == 2)
{
    for(int i=1;i<=b1.m();i++)
    {   
    auto a1=T.cplx(b1(i),St(1));
    auto a2=T.cplx(b1(i),St(3));
    T.set(b1(i),St(1),a2);
    T.set(b1(i),St(3),a1);  
    }
}
if(T.r() == 3)
{   for(int j=1;j <= b2.m();j++)
  {
    for(int i=1;i <= b1.m();i++)
    {   
    auto a1=T.cplx(b1(i),St(1),b2(j));
    auto a2=T.cplx(b1(i),St(3),b2(j));
    T.set(b1(i),St(1),b2(j),a2);
    T.set(b1(i),St(3),b2(j),a1);  
    }
  }
}
   
   psit.setA(k,T);
   //PrintData(psit.A(1));
}

    auto trs = overlapC(psi,psit);
    println("\n TRS check");
    println("\n <psi|psit>",trs);
    println("\n |<psi|psit>|",abs(trs));
/////////////////////////////// end of Time reversal operation that psit is TR of psi
}
MPS
iqmpstomps(IQMPS const& psi, SiteSet const& sites)
{
 auto N = sites.N();
 auto psio = MPS(sites);
 for(int k=1; k <= N; ++k)
 {
	psio.setA(k,toITensor(psi.A(k))); 
 }
 return psio;
}
////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
    {
    //Parse the input file
    if(argc != 2) { printfln("Usage: %s inputfile_exthubbard",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");

    auto Nx = input.getInt("Nx");
    auto Ny = input.getInt("Ny");
    auto hx = input.getReal("hx",0.0);
    Args args;
    args.add("Nx",Nx);
    args.add("Ny",Ny);
    args.add("hx",hx);


    Print(args);

		
//	int Nx = 4;
//    int Ny = 4;
    int N = Nx*Ny;
    auto yperiodic = false;

    //
    // Initialize the site degrees of freedom.
    //
    auto sites = SpinHalf(N);

    //
    // Use the AutoMPO feature to create the 
    // next-neighbor Heisenberg model.
    //
    auto ampo = AutoMPO(sites);
    auto lattice = triangularLattice(Nx,Ny,{"YPeriodic=",yperiodic});
    //square lattice also available:
    //auto lattice = squareLattice(Nx,Ny,{"YPeriodic=",yperiodic});
    for(auto bnd : lattice)
        {
 //       ampo += 0.5,"S+",bnd.s1,"S-",bnd.s2;
//        ampo += 0.5,"S-",bnd.s1,"S+",bnd.s2;
        ampo +=  -1.0,"Sz",bnd.s1,"Sz",bnd.s2;
        }
    for(auto j : range1(N))
        {
		ampo += hx, "Sx", j;
		}	
 
    auto Hq = MPO(ampo);

   

   ITensor Hqa = Hq.A(1);
     for(int j = 2; j <= N; ++j)
        {
         Hqa = Hqa*(Hq.A(j));
        }
println("\n checking point after Hqa");
 //   PrintData(Hqa);


//   PrintData(Hq);
/*
   IQTensor psisa = psis.A(1);
     for(int j = 2; j <= N; ++j)
        {
         psisa = psisa*(psis.A(j));
        }

    PrintData(psisa);
*/
 /*  println("\n For printing Hq");
   PrintData(Hq);
   IQTensor Hqa = Hq.A(1);
     for(int j = 2; j <= N; ++j)
        {
         Hqa = Hqa*(Hq.A(j));
        }

    PrintData(Hqa);
*/
    //
    // overlap(psi,H,psi) = <psi|H|psi>
  /*  auto psis = Singlet(sites);
    IQTensor psisiq = psis.A(1);
     for(int j = 2; j <= N; ++j)
        {
         psisiq = psisiq*(psis.A(j));
        }
    PrintData(psisiq);
    printfln("<psi|S2total|psi> = %.5f", overlap(psis,makeS2(sites),psis) );
*/
    //printfln("<psi|Sx2total|psi> = %.5f", overlap(psis,makeSx2(sites),psis) );
    //printfln("<psi|Sy2total|psi> = %.5f", overlap(psis,makeSy2(sites),psis) );
    //printfln("<psi|Sz2total|psi> = %.5f", overlap(psis,makeSz2(sites),psis) );
 //   printfln("<psis|S2total|psis> = %.5f", overlap(psis,makeS2(sites),psis) );
 //   printfln("Print out S2= %.5f");
 //   PrintData(makeS2(sites));
    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
/*
    auto sweeps = Sweeps(5);
    sweeps.maxm() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    println(sweeps);
*/
    //
    // Begin the DMRG calculation
    //
   // auto energy = dmrg(psi,H,sweeps,"Quiet");

    //
    // Print the final energy reported by DMRG
    //
  //  printfln("\nGround State Energy = %.10f",energy);
    ITensor U,D;

    auto E = diagHermitian(Hqa,U,D);
    println("Eigenvalues: \n");
  //  PrintData(D);
//	print(D);
/*
    auto vinds = std::vector<Index>(2);   
    int iind=0;
    for(auto& I : D.inds())
    {
		vinds[iind]= I;
		iind++;
    }
	auto d1=vinds[0];
	auto d2=vinds[1];
	Print(d1);
	Print(d2);
*/
	for(auto i:range1(int(pow(2,N))))
	{
		printfln("%.10f",D.real(i,i));
	}
	
	std::list<Real> Eenlist = {};
	std::list<Real> betalist = {};
	std::list<Real> templist = {};
	
	for(int iT=0; iT<=200;iT++)
	{
		auto T=iT/100.0+1e-9;
		auto beta=1/T;
		auto Z=0.0;
		auto Een=0.0;
		
		for(auto i:range1(int(pow(2,N))))
		{
			Een=Een+D.real(i,i)*exp(-beta*D.real(i,i));
			Z=Z+exp(-beta*D.real(i,i));
		}
		Een=(Een/Z)/N;
		Eenlist.push_back(Een); 
		betalist.push_back(1.0*beta);
		templist.push_back(1.0/beta); 
	}
	std::ofstream enf("en_beta_temp.dat");
    for(auto n : range(betalist.size()))
        {

			enf << format("%.14f %.14f %.14f\n",betalist.front(),templist.front(),Eenlist.front());
			betalist.pop_front();
			templist.pop_front();
			Eenlist.pop_front();
        }
    enf.close();
//    println("Eigenvectors: \n");
//    PrintData(U);


    return 0;
    }
