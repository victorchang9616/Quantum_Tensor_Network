#include "itensor/all.h"
#include <time.h>
#include <sys/time.h>
#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651
using namespace itensor;
////////////////////////////////////////subroutines//////////////////////////////////////
double
read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}
void
excitedstates(MPS& psiq0,MPS& psiq1,MPS& psiq2,MPS& psiq3,MPO const& Hq,Sweeps const& sweeps)
{
   auto wei1 = 20.0;
   auto wei2 = 20.0;
   auto wei3 = 20.0;
/*
    auto sweeps = Sweeps(5);
    sweeps.maxm() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 3;
    sweeps.noise() = 1E-7,1E-8,0.0;
*/
   auto wfsq = std::vector<MPS>(1);
    wfsq.at(0) = psiq0;
    //auto psi1 = MPS(state);

    //
    // Begin the second state DMRG calculation
    //
  //  auto en1 = dmrg(psi1,H,wfs,sweeps,{"Quiet",true,"Weight",20.0});
    auto en1q = dmrg(psiq1,Hq,wfsq,sweeps,{"Quiet",true,"Weight",wei1});
//    auto psi1 = iqmpstomps(psiq1,sites);
    auto wfsq2 = std::vector<MPS>(2);
    wfsq2.at(0) = psiq0;
    wfsq2.at(1) = psiq1;
    auto otho = overlapC(psiq0,psiq1);
//    println("\n <psiq0|psiq1>",otho);
    println("\n |<psiq0|psiq1>|",abs(otho));
//    auto psiq2 = IQMPS(state);
    auto en2q = dmrg(psiq2,Hq,wfsq2,sweeps,{"Quiet",true,"Weight",wei2});
 //   auto psi2 = iqmpstomps(psiq2,sites);
//    auto psiq3 = IQMPS(state);
    auto wfsq3 = std::vector<MPS>(3);
    wfsq3.at(0) = psiq0;
    wfsq3.at(1) = psiq1;
    wfsq3.at(2) = psiq2;
    auto otho1 = overlapC(psiq1,psiq2);
//    println("\n <psiq2|psiq1>",otho1);
   println("\n |<psiq2|psiq1>|",abs(otho1));
    auto otho2 = overlapC(psiq2,psiq0);
//    println("\n <psiq2|psiq0>",otho2);
   println("\n |<psiq2|psiq0>|",abs(otho2));
    auto en3q = dmrg(psiq3,Hq,wfsq3,sweeps,{"Quiet",true,"Weight",wei3});
    otho1 = overlapC(psiq1,psiq3);
//    println("\n <psiq3|psiq1>",otho1);
   println("\n |<psiq3|psiq1>|",abs(otho1));
    otho2 = overlapC(psiq3,psiq0);
//    println("\n <psiq3|psiq0>",otho2);
   println("\n |<psiq3|psiq0>|",abs(otho2));
   println("\n For psi2/psiq2");
   auto otho3 = overlapC(psiq3,psiq2);
//    println("\n <psiq3|psiq2>",otho3);
   println("\n |<psiq3|psiq2>|",abs(otho3));
    printfln("\nen1 = %.10f",en1q);
//    printfln("\nUsing overlap = %.10fand%.10f", overlap(psi1,H,psi1),overlap(psiq1,Hq,psiq1));
//    println("\nTotal QN of Ground State = ",totalQN(psiq1));

//    println("\n For psi1/psiq1");
    printfln("\nen2 = %.10f",en2q);
//    println("\nTotal QN of Ground State = ",totalQN(psiq2));
;
//   println("\n For psi2/psiq2");
    printfln("\nen3 = %.10f",en3q);
//    println("\nTotal QN of Ground State = ",totalQN(psiq3));

//    auto psi3 = iqmpstomps(psiq3,sites);
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

void
Sicalculate(MPS const& psi, SiteSet const& sites)
{
// Compute magnetization for ITensor
    auto N = sites.N();
    auto psio = psi;
    auto Sxj = 0.0+0.0_i;
    auto Syj = 0.0+0.0_i; 
    auto Szj = 0.0+0.0_i;
    auto Sz2j = 0.0+0.0_i;
    auto Qzzj = 0.0+0.0_i;
    auto Qx2y2j= 0.0+0.0_i;
    auto Qxyj= 0.0+0.0_i;
    auto Qyzj= 0.0+0.0_i;
    auto Qxzj= 0.0+0.0_i;
    auto Sxarray=std::vector<Real>(N);
    auto Syarray=std::vector<Real>(N);
    auto Szarray=std::vector<Real>(N);
    auto Sz2array=std::vector<Real>(N);
    auto Qzzarray=std::vector<Complex>(N);
    auto Qx2y2array=std::vector<Complex>(N);
    auto Qxyarray=std::vector<Complex>(N);
    auto Qyzarray=std::vector<Complex>(N);
    auto Qxzarray=std::vector<Complex>(N);
   FILE *sxf = fopen("sx.txt", "w");
   FILE *syf = fopen("sy.txt", "w");
   FILE *szf = fopen("sz.txt", "w");
   FILE *szf3 = fopen("sz3.txt", "w");
   FILE *ssum = fopen("sxyztot.txt", "w");
    for(int j = 1; j <= N; ++j)
    {
     ITensor Sx = sites.op("Sx",j);
     ITensor Sy = sites.op("Sy",j);
     ITensor Sz = sites.op("Sz",j);
     ITensor Sz2 = sites.op("Sz2",j);
     ITensor Sx2 = sites.op("Sx2",j);
     ITensor Sy2 = sites.op("Sy2",j);
     ITensor SxSy = sites.op("SxSy",j);
     ITensor SySx = sites.op("SySx",j);
     ITensor SxSz = sites.op("SxSz",j);
     ITensor SzSx = sites.op("SzSx",j);
     ITensor SySz = sites.op("SySz",j);
     ITensor SzSy = sites.op("SzSy",j);
	 

    //Make site j the MPS "orthogonality center"
    psio.position(j);
    //Measure magnetization
     auto Sxex = (psio.A(j)* Sx * dag(prime(psio.A(j),Site))).cplx();
     auto Syex = (psio.A(j)* Sy * dag(prime(psio.A(j),Site))).cplx();
     auto Szex = (psio.A(j)* Sz * dag(prime(psio.A(j),Site))).cplx();
     auto Sz2ex = (psio.A(j)* Sz2 * dag(prime(psio.A(j),Site))).cplx();
     auto Sy2ex = (psio.A(j)* Sy2 * dag(prime(psio.A(j),Site))).cplx();
     auto Sx2ex = (psio.A(j)* Sx2 * dag(prime(psio.A(j),Site))).cplx();
     auto SxSyex = (psio.A(j)* SxSy * dag(prime(psio.A(j),Site))).cplx();
     auto SySxex = (psio.A(j)* SySx * dag(prime(psio.A(j),Site))).cplx();
     auto SxSzex = (psio.A(j)* SxSz * dag(prime(psio.A(j),Site))).cplx();
     auto SzSxex = (psio.A(j)* SzSx * dag(prime(psio.A(j),Site))).cplx();
     auto SySzex = (psio.A(j)* SySz * dag(prime(psio.A(j),Site))).cplx();
     auto SzSyex = (psio.A(j)* SzSy * dag(prime(psio.A(j),Site))).cplx();
//     auto SxSyex = (psio.A(j)* Sx*Sy* dag(prime(psio.A(j),Site))).cplx();
     auto Qzz= (3.0*Sz2ex - 2.0)/sqrt(3.0);
     auto Qx2y2= (Sx2ex-Sy2ex);
     auto Qxy= (SxSyex+SySxex);
     auto Qxz= (SxSzex+SzSxex);
     auto Qyz= (SySzex+SzSyex);

//     auto Qxy= (SxSy-Syex);
    Sxarray[j-1]=Sxex.real();
    Syarray[j-1]=Syex.real();
    Szarray[j-1]=Szex.real();
    Sz2array[j-1]=Sz2ex.real();
    Qzzarray[j-1]=Qzz;
    Qx2y2array[j-1]=Qx2y2;
    Qxyarray[j-1]=Qxy;
    Qxzarray[j-1]=Qxz;
    Qyzarray[j-1]=Qyz;
    Sxj = Sxj + Sxex;
    Syj = Syj + Syex;
    Szj = Szj + Szex;
    Sz2j = Sz2j + Sz2ex;
    Qzzj = Qzzj + Qzz;
    Qx2y2j = Qx2y2j + Qx2y2;
    Qxyj = Qxyj + Qxy;
    Qxzj = Qxzj + Qxz;
    Qyzj = Qyzj + Qyz;
   // Syj = Syj + (psi.A(j)* Sy * dag(prime(psi.A(j),Site))).cplx();
    //Szj = Szj + (psi.A(j)* sites.op("Sz",j)* dag(prime(psi.A(j),Site))).cplx();
   //S2j = S2j + (psi.A(j)* S2 * dag(prime(psi.A(j),Site))).cplx();
    //println("Sz_",j," = ",Szj);
//    printf("\nSx%d = %.5f  ",j,Sxex);
//    printf("\nSy%d = %.5f",j,Syex);
//    printf("\nSz%d = %.5f",j,Szex);
//    printf("\nS2%d = %.5f ",j,S2ex);
    }
   println("\n===============");

  fprintf(ssum, "\nrealx=%.10f \nimagx=%.10f", Sxj.real(),Sxj.imag());
  fprintf(ssum, "\nrealy=%.10f \nimagy=%.10f", Syj.real(),Syj.imag());
  fprintf(ssum, "\nrealz=%.10f \nimagz=%.10f", Szj.real(),Szj.imag());
    println("\nSxtotal =  ",Sxj);
    println("\nSytotal =  ",Syj);
    println("\nSztotal =  ",Szj);
    println("\nSz2total =  ",Sz2j);
    println("\nQzztotal =  ",Qzzj);
    println("\n ===============");
    println("\n ===============");
    println("\n local Sx");
     for(int j = 1; j <= N; ++j)
    {
     printf("\n %.10f  ",Sxarray[j-1]);
	fprintf(sxf, "%.10f\n", Sxarray[j-1]);
    }
    println("\n ===============");
    println("\n local Sy");
     for(int j = 1; j <= N; ++j)
    {
     printf("\n %.10f  ",Syarray[j-1]);
	fprintf(syf, "%.10f\n", Syarray[j-1]);
    }
    println("\n ===============");
    println("\n local Sz");
     for(int j = 1; j <= N; ++j)
    {
     printf("\n %.10f  ",Szarray[j-1]);
	fprintf(szf, "%.10f\n", Szarray[j-1]);
    }
    println("\n ===============");
    println("\n local Sz3");
     for(int j=1;j<= N;j=j+3)
    { 
     printf("\n %.10f  ",Szarray[j-1]);
        fprintf(szf3, "%.10f\n", Szarray[j-1]);
    }
   println("\n ==============="); 
    println("\n local Sz*Sz");
     for(int j = 1; j <= N; ++j)
    {
     printf("\n %.10f  ",Sz2array[j-1]);
//	fprintf(szf, "%.10f\n", Szarray[j-1]);
    }
   println("\n ==============="); 
    println("\n local Qzz");
     for(int j = 1; j <= N; ++j)
    {
     printf("\n %.10f  ",Qzzarray[j-1].real());
//	fprintf(szf, "%.10f\n", Szarray[j-1]);
    }
   println("\n ==============="); 
    println("\n local Qzz3");
     for(int j=1;j<= N;j=j+3)
    { 
     printf("\n %.10f  ",Qzzarray[j-1].real());
  //      fprintf(szf3, "%.10f\n", Szarray[j-1]);
    }
   println("\n ==============="); 
    println("\n local Qx2y2");
     for(int j=1;j<= N;j=j+1)
    { 
     printf("\n %.10f ",Qx2y2array[j-1].real());
  //      fprintf(szf3, "%.10f\n", Szarray[j-1]);
    }
   println("\n ==============="); 
    println("\n local Qx2y23");
	     for(int j=1;j<= N;j=j+3)
    { 
     printf("\n %.10f",Qx2y2array[j-1].real());
  //      fprintf(szf3, "%.10f\n", Szarray[j-1]);
    }
    println("\n local Qxy");
	     for(int j=1;j<= N;j=j+1)
    { 
     printf("\n%.10f",Qxyarray[j-1].real());
  //      fprintf(szf3, "%.10f\n", Szarray[j-1]);
    }
	    println("\n local Qxy3");
	     for(int j=1;j<= N;j=j+3)
    { 
     printf("\n%.10f",Qxyarray[j-1].real());
  //      fprintf(szf3, "%.10f\n", Szarray[j-1]);
    }
    println("\n local Qxz");
	     for(int j=1;j<= N;j=j+1)
    { 
     printf("\n%.10f",Qxzarray[j-1].real());
  //      fprintf(szf3, "%.10f\n", Szarray[j-1]);
    }
	    println("\n local Qxz3");
	     for(int j=1;j<= N;j=j+3)
    { 
     printf("\n%.10f",Qxzarray[j-1].real());
  //      fprintf(szf3, "%.10f\n", Szarray[j-1]);
    }
    println("\n local Qxy");
	     for(int j=1;j<= N;j=j+1)
    { 
     printf("\n%.10f",Qxyarray[j-1].real());
  //      fprintf(szf3, "%.10f\n", Szarray[j-1]);
    }
	    println("\n local Qxy3");
	     for(int j=1;j<= N;j=j+3)
    { 
     printf("\n%.10f",Qxyarray[j-1].real());
  //      fprintf(szf3, "%.10f\n", Szarray[j-1]);
    }
    println("\n local Qyz");
	     for(int j=1;j<= N;j=j+1)
    { 
     printf("\n%.10f",Qyzarray[j-1].real());
  //      fprintf(szf3, "%.10f\n", Szarray[j-1]);
    }
	    println("\n local Qyz3");
	     for(int j=1;j<= N;j=j+3)
    { 
     printf("\n%.10f",Qyzarray[j-1].real());
  //      fprintf(szf3, "%.10f\n", Szarray[j-1]);
    }
   println("\n ==============="); 
   fclose(sxf);
   fclose(syf);
   fclose(szf);
   fclose(ssum);
   fclose(szf3);

}
void
trimerbond(double K11,double K12,double K21,double K22,MPS const& psi, SiteSet const& sites)
{
// Compute magnetization for ITensor
    auto N = sites.N();
    auto psio = psi;
    auto Hnn=std::vector<Complex>(N);
    auto Hnnn=std::vector<Complex>(N);
	FILE *Hnnf = fopen("Hnn.txt", "w");
	FILE *Hnnnf = fopen("Hnnn.txt", "w");
/*
    auto Sxj = 0.0+0.0_i;
    auto Syj = 0.0+0.0_i; 
    auto Szj = 0.0+0.0_i;
    auto S2j = 0.0+0.0_i;
*/

for(int i1 = 1; i1 <= (N-1); ++i1)
{
        //Print(psi.A(i1-1));
        //Print(psi.A(i1));
        //Print(psi.A(i1+1));
        psio.position(i1);
       
        //ITensor Sx1 = 0.5*(sites.op("Sp",i1) + sites.op("Sm",i1));
        //ITensor Sy1 = (1/(2_i))*(sites.op("Sp",i1) - sites.op("Sm",i1));
        //ITensor Sz1 = sites.op("Sz",i1);

        ITensor Sp1 = sites.op("Sp",i1);
        ITensor Sm1 = sites.op("Sm",i1);
        ITensor Sz1 = sites.op("Sz",i1);
        ITensor Sp2 = sites.op("Sp",i1+1);
        ITensor Sm2 = sites.op("Sm",i1+1);
        ITensor Sz2 = sites.op("Sz",i1+1);
/*
        ITensor Sx3 = sites.op("Sp",i1+2);
        ITensor Sy3 = sites.op("Sm",i1+2);
        ITensor Sz3 = sites.op("Sz",i1+2);
*/

        auto B1 = (K11/4.0)*Sp1*Sp1*Sm2*Sm2 + (K11/4.0)*Sm1*Sm1*Sp2*Sp2 + (K11/4.0)*Sp1*Sm1*Sm2*Sp2 + (K11/4.0)*Sm1*Sp1*Sp2*Sm2 
		+ (K11/1.0)*Sz1*Sz1*Sz2*Sz2 
		+ (K11/2.0)*Sp1*Sz1*Sm2*Sz2 + (K11/2.0)*Sm1*Sz1*Sp2*Sz2 + (K11/2.0)*Sz1*Sp1*Sz2*Sm2	+ (K11/2.0)*Sz1*Sm1*Sz2*Sp2
		+ (0.5)*K12*Sp1*Sm2 + (0.5)*K12*Sm1*Sp2 + K12*Sz1*Sz2;
        auto wf1 = psio.A(i1)*psio.A(i1+1);
        // compute <wf1|B1|wf1>
        Hnn[i1-1]=(dag(prime(wf1,Site)) * B1 * wf1).cplx();
//        Kaiq += kai[i1-1];
}
for(int i1 = 1; i1 <= (N-2); ++i1)
{
        //Print(psi.A(i1-1));
        //Print(psi.A(i1));
        //Print(psi.A(i1+1));
        psio.position(i1);
       
        //ITensor Sx1 = 0.5*(sites.op("Sp",i1) + sites.op("Sm",i1));
        //ITensor Sy1 = (1/(2_i))*(sites.op("Sp",i1) - sites.op("Sm",i1));
        //ITensor Sz1 = sites.op("Sz",i1);

        ITensor Sp1 = sites.op("Sp",i1);
        ITensor Sm1 = sites.op("Sm",i1);
        ITensor Sz1 = sites.op("Sz",i1);
        ITensor Sp2 = sites.op("Sp",i1+2);
        ITensor Sm2 = sites.op("Sm",i1+2);
        ITensor Sz2 = sites.op("Sz",i1+2);
        auto B1 = (K21/4.0)*Sp1*Sp1*Sm2*Sm2 + (K21/4.0)*Sm1*Sm1*Sp2*Sp2 + (K21/4.0)*Sp1*Sm1*Sm2*Sp2 + (K21/4.0)*Sm1*Sp1*Sp2*Sm2 
		+ (K21/1.0)*Sz1*Sz1*Sz2*Sz2 
		+ (K21/2.0)*Sp1*Sz1*Sm2*Sz2 + (K21/2.0)*Sm1*Sz1*Sp2*Sz2 + (K21/2.0)*Sz1*Sp1*Sz2*Sm2	+ (K21/2.0)*Sz1*Sm1*Sz2*Sp2
		+ (0.5)*K22*Sp1*Sm2 + (0.5)*K22*Sm1*Sp2 + K22*Sz1*Sz2;
        auto wf1 = psio.A(i1)*psio.A(i1+2);
        // compute <wf1|B1|wf1>
        Hnnn[i1-1]=(dag(prime(wf1,Site)) * B1 * wf1).cplx();
//        Kaiq += kai[i1-1];
}
//  printfln("\nKaiq = %.10f ",Kaiq);
//  fprintf(kaisumf, "\nreal=%.10f \nimag=%.10f", Kaiq.real(),Kaiq.imag());
  println("\n ========Hnn=======");
    for(int j = 1; j <= N; ++j)
    {
     printf("\n %.10f  ",Hnn[j-1].real());
	 fprintf(Hnnf, "%.10f\n", Hnn[j-1].real());
    }
   println("\n =======Hnnn=====");
	for(int j = 1; j <= N; ++j)
    {
     printf("\n %.10f  ",Hnnn[j-1].real());
	 fprintf(Hnnnf, "%.10f\n", Hnnn[j-1].real());
    }
  println("\n ===============");
  fclose(Hnnf);
  fclose(Hnnnf);

}
void
Kaicalculate(MPS const& psi,  SiteSet const& sites)
{
//  Calculate order of pararmeter Kai = Sumover i { S_i dot S_i+1 cross S_i+2 (start)}
   auto N = sites.N();
   auto psio = psi;
   auto Kaiq= 0.0+0.0_i;
	auto kai=std::vector<Complex>(N);
	FILE *kaif = fopen("kaiF.txt", "w");
	FILE *kaisumf = fopen("kaiFsum.txt", "w");

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
        auto wf1 = psio.A(i1)*psio.A(i1+1)*psio.A(i1+2);
        // compute <wf1|B1|wf1>
        kai[i1-1]=(dag(prime(wf1,Site)) * B1 * wf1).cplx();
        Kaiq += kai[i1-1];
}
  printfln("\nKaiq = %.10f ",Kaiq);
  fprintf(kaisumf, "\nreal=%.10f \nimag=%.10f", Kaiq.real(),Kaiq.imag());
  println("\n ========Fchiral=======");
    for(int j = 1; j <= N; ++j)
    {
     printf("\n %.10f  ",kai[j-1].real());
	 fprintf(kaif, "%.10f\n", kai[j-1].real());
    }
  println("\n ===============");
  fclose(kaif);
  fclose(kaisumf);
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
////////////////////////////////////////subroutines//////////////////////////////////////
int main(int argc, char* argv[])
{
    if(argc != 2) 
   { 
   //reminds us to give an input file if we forget
   printfln("Usage: %s inputfile",argv[0]); 
   return 0; 
   }
   auto input = InputGroup(argv[1],"input");
  auto N = input.getInt("N"); //number of sites
  auto Nex = input.getInt("Nex"); //number of sites
  auto h = input.getReal("h",0.0); //coupling
   auto qzz = input.getReal("qzz",0.0); //coupling
  auto qx2y2 = input.getReal("qx2y2",0.0); //coupling
  auto qxy = input.getReal("qxy",0.0); //coupling
  auto qyz = input.getReal("qyz",0.0); //coupling
  auto qxz = input.getReal("qxz",0.0); //coupling
  auto g = input.getReal("g",0.0); //coupling
  auto wei1 = input.getReal("wei1",20.0); //coupling
  auto lamda = input.getReal("lamda",0.0); //coupling
  auto f = input.getReal("f",0.0); //coupling
  auto Jp = input.getReal("Jp",-1.0); //coupling
  auto Kone = input.getReal("Kone",1.0); //coupling
   auto Ktwo = input.getReal("Ktwo",1.0); //coupling
//   auto K11 = input.getReal("K11",1.0); //coupling
//   auto K12 = input.getReal("K12",1.0); //coupling
   auto angleJp = input.getReal("angleJp",90.0); //coupling
   auto angleK = input.getReal("angleK",90.0); //coupling
 
//    int N = 63;

    //
    // Initialize the site degrees of freedom.
    //
    //auto sites = SpinHalf(N); //make a chain of N spin 1/2's
 //   SpinOne sites;
 //   readFromFile("sites_318",sites);
    auto sites = SpinOne(N); //make a chain of N spin 1's
// parameter tuning 
//    double angle = 20.0;
 //   double lamda = 0.3;
 //   double J = 1.0;
// parameter tuning
//    double K = 1.0;
    double thetaK = angleK*PI/180.0;
    double K11 = Kone*cos(thetaK);
    double K12 = Kone*sin(thetaK);
    double K21 = Ktwo*cos(thetaK);
    double K22 = Ktwo*sin(thetaK);
    double thetaJp = angleJp*PI/180.0;
    double Jp1 = Jp*cos(thetaJp);
    double Jp2 = Jp*sin(thetaJp);
 //   double Jp = J*sin(PI/8.0);
 //   double Jp = -1.0;

 //   printfln("Kone = %.10f \nKtwo = %.10f\n Jp = %.10f\nh = %.10f\nN= %i\n ",Kone,Ktwo,Jp,h,N);
 //   printfln("lamda=%.10f\n",lamda);
    //
    // Use the AutoMPO feature to create the 
    // next-neighbor Heisenberg model.
    //
    // Here we convert the AutoMPO information
    // into an IQMPO, a matrix-product operator
    // which automatically tracks quantum
    // number information.
    //printmatrixstateq(IQMPS const& psiq, int N)

    auto ampo = AutoMPO(sites);// for H1D
   // auto ampo1 = AutoMPO(sites);//for Hp
    auto ampoKone = AutoMPO(sites);
    auto ampoKtwo = AutoMPO(sites);
    auto ampoJp = AutoMPO(sites);
    auto ampolamdap = AutoMPO(sites);
    auto ampolamdam = AutoMPO(sites);
    auto amposz = AutoMPO(sites);
    auto ampof = AutoMPO(sites);
   // auto ampo1 = AutoMPO(sites);//for Hp


     for(int j = 1; j <= (N-1); ++j)
        {
 
         ampo += K11/4.0,"Sp*Sp",j,"Sm*Sm",j+1;
         ampo += K11/4.0,"Sm*Sm",j,"Sp*Sp",j+1;
         ampo += K11/4.0,"Sp*Sm",j,"Sm*Sp",j+1;
         ampo += K11/4.0,"Sm*Sp",j,"Sp*Sm",j+1;
         ampo += K11/1.0,"Sz*Sz",j,"Sz*Sz",j+1;
         ampo += K11/2.0,"Sp*Sz",j,"Sm*Sz",j+1;
         ampo += K11/2.0,"Sm*Sz",j,"Sp*Sz",j+1;
         ampo += K11/2.0,"Sz*Sp",j,"Sz*Sm",j+1;
         ampo += K11/2.0,"Sz*Sm",j,"Sz*Sp",j+1;
 
         ampo += (0.5)*K12,"S+",j,"S-",j+1;
         ampo += (0.5)*K12,"S-",j,"S+",j+1;
         ampo += (1.0)*K12,"Sz",j,"Sz",j+1;
         }



     for(int j = 1; j <= (N-2); ++j)
        {
         ampo += K21/4.0,"Sp*Sp",j,"Sm*Sm",j+2;
         ampo += K21/4.0,"Sm*Sm",j,"Sp*Sp",j+2;
         ampo += K21/4.0,"Sp*Sm",j,"Sm*Sp",j+2;
         ampo += K21/4.0,"Sm*Sp",j,"Sp*Sm",j+2;
         ampo += K21/1.0,"Sz*Sz",j,"Sz*Sz",j+2;
         ampo += K21/2.0,"Sp*Sz",j,"Sm*Sz",j+2;
         ampo += K21/2.0,"Sm*Sz",j,"Sp*Sz",j+2;
         ampo += K21/2.0,"Sz*Sp",j,"Sz*Sm",j+2;
         ampo += K21/2.0,"Sz*Sm",j,"Sz*Sp",j+2;
         ampo += (0.5)*K22,"S+",j,"S-",j+2;
         ampo += (0.5)*K22,"S-",j,"S+",j+2;
         ampo += (1.0)*K22,"Sz",j,"Sz",j+2;
        }
      for(int j = 1; j <= (N-3); ++j)
        {
       ampo += -Jp1/4.0,"Sp*Sp",j,"Sm*Sm",j+3;
         ampo += -Jp1/4.0,"Sm*Sm",j,"Sp*Sp",j+3;
         ampo += -Jp1/4.0,"Sp*Sm",j,"Sm*Sp",j+3;
         ampo += -Jp1/4.0,"Sm*Sp",j,"Sp*Sm",j+3;
         ampo += -Jp1/1.0,"Sz*Sz",j,"Sz*Sz",j+3;
         ampo += -Jp1/2.0,"Sp*Sz",j,"Sm*Sz",j+3;
         ampo += -Jp1/2.0,"Sm*Sz",j,"Sp*Sz",j+3;
         ampo += -Jp1/2.0,"Sz*Sp",j,"Sz*Sm",j+3;
         ampo += -Jp1/2.0,"Sz*Sm",j,"Sz*Sp",j+3;       
         ampo += (-0.5)*Jp2,"S+",j,"S-",j+3;
         ampo += (-0.5)*Jp2,"S-",j,"S+",j+3;
         ampo += (-1.0)*Jp2,"Sz",j,"Sz",j+3;
         }
     

// for ferrDSSCO Hp > 0
/*
      for(int j = 1; j <= (N-2); ++j)
        {
         ampo += (0.5_i)*lamda,"Sz",j,"Sp",j+1,"Sm",j+2;
         ampo += -(0.5_i)*lamda,"Sz",j,"Sm",j+1,"Sp",j+2;
         ampo += -(0.5_i)*lamda,"Sp",j,"Sz",j+1,"Sm",j+2;
         ampo += (0.5_i)*lamda,"Sm",j,"Sz",j+1,"Sp",j+2;
         ampo += (0.5_i)*lamda,"Sp",j,"Sm",j+1,"Sz",j+2;
         ampo += -(0.5_i)*lamda,"Sm",j,"Sp",j+1,"Sz",j+2;

        }
*/
// for ferrDSSCO Hp > 0

// for AFDSSCO Hp > 0 if i is even



// for AFDSSCO Hp > 0
// for AFDSSCO Hp < 0 if i is odd

     for(int j = 1; j <= 1; ++j)
        {
 
         ampo += h,"Sz",j;

         }
       for(int j = 1; j <= (N-2); j=j+1)
        {
         ampo += -(0.5_i)*lamda,"Sz",j,"Sp",j+1,"Sm",j+2;
         ampo += (0.5_i)*lamda,"Sz",j,"Sm",j+1,"Sp",j+2;
         ampo += (0.5_i)*lamda,"Sp",j,"Sz",j+1,"Sm",j+2;
         ampo += -(0.5_i)*lamda,"Sm",j,"Sz",j+1,"Sp",j+2;
         ampo += -(0.5_i)*lamda,"Sp",j,"Sm",j+1,"Sz",j+2;
         ampo += (0.5_i)*lamda,"Sm",j,"Sp",j+1,"Sz",j+2;

        }
       for(int j = 1; j <= (1); j=j+1)
        {
         ampo += (2.0)*qzz,"Sz*Sz",j;
	}
       for(int j = 1; j <= (1); j=j+1)
        {
         ampo += (1.0/2.0)*qx2y2,"Sp*Sp",j;
         ampo += (1.0/2.0)*qx2y2,"Sm*Sm",j;
		}
       for(int j = 1; j <= (1); j=j+1)
        {
         ampo += (1.0/2.0_i)*qxy,"Sp*Sp",j;
         ampo += (-1.0/2.0_i)*qxy,"Sm*Sm",j;
		}
       for(int j = 1; j <= (1); j=j+1)
        {
         ampo += (1.0/2.0_i)*qyz,"Sp*Sz",j;
         ampo += (-1.0/2.0_i)*qyz,"Sm*Sz",j;
         ampo += (1.0/2.0_i)*qyz,"Sz*Sp",j;
         ampo += (-1.0/2.0_i)*qyz,"Sz*Sm",j;
		}
       for(int j = 1; j <= (1); j=j+1)
        {
         ampo += (1.0/2.0)*qxz,"Sp*Sz",j;
         ampo += (1.0/2.0)*qxz,"Sm*Sz",j;
         ampo += (1.0/2.0)*qxz,"Sz*Sp",j;
         ampo += (1.0/2.0)*qxz,"Sz*Sm",j;
		}

       for(int j = 1; j <= (1); j=j+1)
        {
         ampo += (2.0)*qzz,"Sz*Sz",j;
	}
// for AFDSSCO Hp < 0

 for(int j = 1; j <= (N-1); ++j)
        {
 
         ampoKone += K11/4.0,"Sp*Sp",j,"Sm*Sm",j+1;
         ampoKone += K11/4.0,"Sm*Sm",j,"Sp*Sp",j+1;
         ampoKone += K11/4.0,"Sp*Sm",j,"Sm*Sp",j+1;
         ampoKone += K11/4.0,"Sm*Sp",j,"Sp*Sm",j+1;
         ampoKone += K11/1.0,"Sz*Sz",j,"Sz*Sz",j+1;
         ampoKone += K11/2.0,"Sp*Sz",j,"Sm*Sz",j+1;
         ampoKone += K11/2.0,"Sm*Sz",j,"Sp*Sz",j+1;
         ampoKone += K11/2.0,"Sz*Sp",j,"Sz*Sm",j+1;
         ampoKone += K11/2.0,"Sz*Sm",j,"Sz*Sp",j+1;
 
         ampoKone += (0.5)*K12,"S+",j,"S-",j+1;
         ampoKone += (0.5)*K12,"S-",j,"S+",j+1;
         ampoKone += (1.0)*K12,"Sz",j,"Sz",j+1;
         }



     for(int j = 1; j <= (N-2); ++j)
        {
         ampoKtwo += K21/4.0,"Sp*Sp",j,"Sm*Sm",j+2;
         ampoKtwo += K21/4.0,"Sm*Sm",j,"Sp*Sp",j+2;
         ampoKtwo += K21/4.0,"Sp*Sm",j,"Sm*Sp",j+2;
         ampoKtwo += K21/4.0,"Sm*Sp",j,"Sp*Sm",j+2;
         ampoKtwo += K21/1.0,"Sz*Sz",j,"Sz*Sz",j+2;
         ampoKtwo += K21/2.0,"Sp*Sz",j,"Sm*Sz",j+2;
         ampoKtwo += K21/2.0,"Sm*Sz",j,"Sp*Sz",j+2;
         ampoKtwo += K21/2.0,"Sz*Sp",j,"Sz*Sm",j+2;
         ampoKtwo += K21/2.0,"Sz*Sm",j,"Sz*Sp",j+2;
         ampoKtwo += (0.5)*K22,"S+",j,"S-",j+2;
         ampoKtwo += (0.5)*K22,"S-",j,"S+",j+2;
         ampoKtwo += (1.0)*K22,"Sz",j,"Sz",j+2;
 //        ampoKtwo += (0.5)*Ktwo,"S+",j,"S-",j+2;
 //        ampoKtwo += (0.5)*Ktwo,"S-",j,"S+",j+2;
 //        ampoKtwo += (1.0)*Ktwo,"Sz",j,"Sz",j+2;
        }
      for(int j = 1; j <= (N-3); ++j)
        {
         ampoJp += -Jp1/4.0,"Sp*Sp",j,"Sm*Sm",j+3;
         ampoJp += -Jp1/4.0,"Sm*Sm",j,"Sp*Sp",j+3;
         ampoJp += -Jp1/4.0,"Sp*Sm",j,"Sm*Sp",j+3;
         ampoJp += -Jp1/4.0,"Sm*Sp",j,"Sp*Sm",j+3;
         ampoJp += -Jp1/1.0,"Sz*Sz",j,"Sz*Sz",j+3;
         ampoJp += -Jp1/2.0,"Sp*Sz",j,"Sm*Sz",j+3;
         ampoJp += -Jp1/2.0,"Sm*Sz",j,"Sp*Sz",j+3;
         ampoJp += -Jp1/2.0,"Sz*Sp",j,"Sz*Sm",j+3;
         ampoJp += -Jp1/2.0,"Sz*Sm",j,"Sz*Sp",j+3;       
         ampoJp += (-0.5)*Jp2,"S+",j,"S-",j+3;
         ampoJp += (-0.5)*Jp2,"S-",j,"S+",j+3;
         ampoJp += (-1.0)*Jp2,"Sz",j,"Sz",j+3;
         }
      for(int j = 1; j <= 1; ++j)
        {
 
         amposz += h,"Sz",j;

         }
     for(int j = 1; j <= N; ++j)
        {
 
         ampof += f,"Sz",j,"Sz",j;

         }
//    Hp for FerroDSSCO Hp >0 for all i
/*     
      for(int j = 1; j <= (N-2); ++j)
        {
         ampolamda += (0.5_i)*lamda,"Sz",j,"Sp",j+1,"Sm",j+2;
         ampolamda += -(0.5_i)*lamda,"Sz",j,"Sm",j+1,"Sp",j+2;
         ampolamda += -(0.5_i)*lamda,"Sp",j,"Sz",j+1,"Sm",j+2;
         ampolamda += (0.5_i)*lamda,"Sm",j,"Sz",j+1,"Sp",j+2;
         ampolamda += (0.5_i)*lamda,"Sp",j,"Sm",j+1,"Sz",j+2;
         ampolamda += -(0.5_i)*lamda,"Sm",j,"Sp",j+1,"Sz",j+2;

        }
*/
//    Hp for FerroDSSCO Hp >0 if i is even
     
 /*     for(int j = (N-2); j <= (N-2); j=j+2)
        {
	    if(j%2 == 0)
		{
         ampolamdap += (0.5_i)*lamda,"Sz",j,"Sp",j+1,"Sm",j+2;
         ampolamdap += -(0.5_i)*lamda,"Sz",j,"Sm",j+1,"Sp",j+2;
         ampolamdap += -(0.5_i)*lamda,"Sp",j,"Sz",j+1,"Sm",j+2;
         ampolamdap += (0.5_i)*lamda,"Sm",j,"Sz",j+1,"Sp",j+2;
         ampolamdap += (0.5_i)*lamda,"Sp",j,"Sm",j+1,"Sz",j+2;
         ampolamdap += -(0.5_i)*lamda,"Sm",j,"Sp",j+1,"Sz",j+2;		
		}			
        else
		{
         ampolamdam += -(0.5_i)*lamda,"Sz",j,"Sp",j+1,"Sm",j+2;
         ampolamdam += (0.5_i)*lamda,"Sz",j,"Sm",j+1,"Sp",j+2;
         ampolamdam += (0.5_i)*lamda,"Sp",j,"Sz",j+1,"Sm",j+2;
         ampolamdam += -(0.5_i)*lamda,"Sm",j,"Sz",j+1,"Sp",j+2;
         ampolamdam += -(0.5_i)*lamda,"Sp",j,"Sm",j+1,"Sz",j+2;
         ampolamdam += (0.5_i)*lamda,"Sm",j,"Sp",j+1,"Sz",j+2;			
		}
          


        }
*/
//    Hp for FerroDSSCO Hp <0 if i is odd
     
      for(int j = 1; j <= N-2; j=j+1)
        {
         ampolamdam += -(0.5_i)*lamda,"Sz",j,"Sp",j+1,"Sm",j+2;
         ampolamdam += (0.5_i)*lamda,"Sz",j,"Sm",j+1,"Sp",j+2;
         ampolamdam += (0.5_i)*lamda,"Sp",j,"Sz",j+1,"Sm",j+2;
         ampolamdam += -(0.5_i)*lamda,"Sm",j,"Sz",j+1,"Sp",j+2;
         ampolamdam += -(0.5_i)*lamda,"Sp",j,"Sm",j+1,"Sz",j+2;
         ampolamdam += (0.5_i)*lamda,"Sm",j,"Sp",j+1,"Sz",j+2;

        }

    auto Hsum = std::vector<MPO>(5);
    Hsum.at(0) = MPO(ampoKone);
    Hsum.at(1) = MPO(ampoKtwo);
    Hsum.at(2) = MPO(ampoJp);
    Hsum.at(3) = MPO(amposz);
    Hsum.at(4) = MPO(ampolamdam);       
    auto Hq = MPO(ampo);

    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
   auto state = InitState(sites);
/*
for(int i = 1; i <= N; ++i) 
        {
         state.set(i,"Z0");
 
        }
*/
for(int i=1;i<=N;i++)
{
if(i%2 ==  1)
   state.set(i,"Up");
else
   state.set(i,"Dn");
}
/*
    for(int i = 1; i <= N; ++i) 
      {
       if(i%3 == 0)
            state.set(i,"Up");
        }
*/
  //  state.set(1,"Z0");
  //  state.set(2,"Z0");
  //  state.set(3,"Z0");
  //  state.set(4,"Z0");
   // state.set(5,"Z0");
     auto psi0 = MPS(state);
//	readFromFile("psi7_318",psi0);
     auto psi1 = MPS(state);
     auto psi2 = MPS(state);
     auto psi3 = MPS(state);
     auto psi4 = MPS(state);
     auto psi5 = MPS(state);
     auto psi6 = MPS(state);
     auto psi7 = MPS(state);
//     auto psi2 = MPS(state);
//     auto psi3 = MPS(state);
    //PrintData(state);
    //PrintData(psi);
    //
    // overlap calculates matrix elements of MPO's with respect to MPS's
    // overlap(psi,H,psi) = <psi|H|psi>
    //
 //   printfln("Initial energy = %.5f", overlap(psi,H,psi) );

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
    auto nsweeps5 = input.getReal("nsweeps5",5);
    auto sw_table5 = InputGroup(input,"sweeps5");
 //   auto quiet = input.getYesNo("quiet",true);
    auto sweeps5 = Sweeps(nsweeps5,sw_table5);
	
 //   auto quiet = input.getYesNo("quiet",true);
	
 //   auto quiet = input.getYesNo("quiet",true);

 double simulation_time = read_timer( );
//    auto sweeps = Sweeps(nsweeps);
//    /*
//        sweeps.maxm() = 10,20,100,200,200,400,400,800;
//            sweeps.cutoff() = 1E-10;
//                sweeps.niter() = 2;
//                    sweeps.noise() = 1E-7,1E-8,0.0,0.0,0.0;
//                    */
//                        println(sweeps);/
    // Begin the first state DMRG calculation
    //
 //   auto en0 = dmrg(psi0,H,sweeps,"Quiet");
    auto en0 = dmrg(psi0,Hq,sweeps5,{"Quiet=",true,"WriteM=",400});
//	println("\nTotal QN of Ground State = ",totalQN(psiq0));
 //   auto psi0 = iqmpstomps(psiq0,sites);
    // wfs MPS/IQMPS to penalize dmrg energy
  //  auto wfs = std::vector<MPS>(1);
   // wfs.at(0) = psi0;
    writeToFile(format("sites_%d",N),sites); //file name will be sites_100
//    writeToFile(format("psi5_%d",N),psiq0);     //file name will be psi_100
	
//	println("\nTotal QN of Ground State = ",totalQN(psiq0));
//    writeToFile(format("psi6_%d",N),psiq0);     //file name will be psi_100
	
//	println("\nTotal QN of Ground State = ",totalQN(psiq0));
    writeToFile(format("psi0_%d",N),psi0);     //file name will be psi_100
//    printfln("\n ============= %i========",(N-m));
//    printfln("\nEo= %.10f",en0q);
 //   printfln("\nUsing overlap = %.10fand%.10f", overlap(psi0,H,psi0),overlap(psiq0,Hq,psiq0));

//    println("\nTotal QN of Ground State = ",totalQN(psiq0));
    auto S2MPO = makeS2mpo(sites);
    auto S20 = overlap(psi0,S2MPO,psi0);
    auto S0=(-1.0+sqrt(1+4*S20))/2.0;
    printfln("\nen0 = %.10f, S2 = %.10f",en0,S20);
	printfln("\nS = %.10f",S0);
    Kaicalculate(psi0,sites);
//    trimerbond(K11,K12,K21,K22,psi0,sites);
    Sicalculate(psi0,sites);
    println(sweeps5);
   auto wfsq = std::vector<MPS>(1);
    wfsq.at(0) = psi0;
    //auto psi1 = MPS(state);

    //
    // Begin the second state DMRG calculation
    //
  //  auto en1 = dmrg(psi1,H,wfs,sweeps,{"Quiet",true,"Weight",20.0});
    auto en1 = dmrg(psi1,Hq,wfsq,sweeps5,{"Quiet",true,"Weight",wei1});
    writeToFile(format("psi1_%d",N),psi1);
//    auto psi1 = iqmpstomps(psiq1,sites);
    auto wfsq2 = std::vector<MPS>(2);
    wfsq2.at(0) = psi0;
    wfsq2.at(1) = psi1;
    auto S21 = overlap(psi0,S2MPO,psi0);
    auto S1=(-1.0+sqrt(1+4*S21))/2.0;
    printfln("\nen1 = %.10f, S2 = %.10f",en1,S21);
	printfln("\nS = %.10f",S1);
    Kaicalculate(psi1,sites);
//    trimerbond(K11,K12,K21,K22,psi0,sites);
    Sicalculate(psi1,sites);
    auto otho = overlapC(psi0,psi1);
    println("\n |<psiq0|psiq1>|",abs(otho));
       auto psid = std::vector<MPS>(50);
    psid.at(0) = psi0;
    psid.at(1) = psi1; 
    for(int d=2;d<=49;d++) psid.at(d) = MPS(state);

   for(int d=2;d<=Nex;d++)
   {

    auto wfsqd = std::vector<MPS>(d);
    for(int i=d;i>0;i--) wfsqd.at(d-i)= psid.at(d-i);    
    auto en2 = dmrg(psid.at(d),Hq,wfsqd,sweeps5,{"Quiet",true,"Weight",wei1,"WriteM=",400});     
    writeToFile(format("psi%d_%d",d,N),psid.at(d));
//    auto psi1 = iqmpstomps(psiq1,sites);


    auto S22 = overlap(psid.at(d),S2MPO,psid.at(d));
    auto S2=(-1.0+sqrt(1+4*S22))/2.0;
    printfln("\nen%d = %.10f, S2 = %.10f",d,en2,S22);
    printfln("\nS = %.10f",S2);
    Kaicalculate(psid.at(d),sites);
//    trimerbond(K11,K12,K21,K22,psi0,sites);
    Sicalculate(psid.at(d),sites);
    
    for(int i=(d-1);i>=0;i--)
    {    
    auto otho20 = overlapC(psid.at(d),psid.at(i));//check the orthognality
    printf("\n |<psiq%d|psiq%d>|=%.10f",i,d,abs(otho20));
    }
   }
 //   excitedstates(psi0,psi1,psi2,psi3,Hq,sweeps);
 //   printfln("\nUsing overlap = %.10fand%.10f", overlap(psi0,H,psi0),overlap(psiq0,Hq,psiq0));
 //   println("\nTotal QN of Ground State = ",totalQN(psiq0));

    //timereversal(psi0,sites);
   // Sicalculate(psiq0,sites);

    // always save sitefile and mpsfile
//    writeToFile("sitefile",sites);
//    writeToFile("mpsfile",psi0);
/*
    // check TRS of psi////////////////////////////////////
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
    println("\n <psi|psit>",trs);
    println("\n |<psi|psit>|",abs(trs));
/////////////////////////////// end of Time reversal operation that psit is TR of psi
*/
/*
// Compute magnetization (start)
    auto  Sxj = 0.0+0.0_i;
    auto Syj = 0.0+0.0_i; 
    auto Szj = 0.0+0.0_i;
    for(int j = 1; j <= N; ++j)
    {
     ITensor Sx = sites.op("Sx",j);
     ITensor Sy = sites.op("Sy",j);
    //Make site j the MPS "orthogonality center"
    psi.position(j);
    //Measure magnetization
    //Sxj = Sxj + (psi.A(j)* sites.op("Sx",j)* dag(prime(psi.A(j),Site))).real();
    Sxj = Sxj + (psi.A(j)* Sx * dag(prime(psi.A(j),Site))).cplx();
    Syj = Syj + (psi.A(j)* Sy * dag(prime(psi.A(j),Site))).cplx();
    //Syj = Syj + (psi.A(j)* sites.op("Sy",j)* dag(prime(psi.A(j),Site))).cplx();
    Szj = Szj + (psi.A(j)* sites.op("Sz",j)* dag(prime(psi.A(j),Site))).cplx();
    //println("Sz_",j," = ",Szj);
    }
    println("\nSxtotal =  ",Sxj);
    println("\nSytotal =  ",Syj);
    println("\nSztotal =  ",Szj);
    // (end)
*/
//  Calculate order of pararmeter Kai = Sumover i { S_i dot S_i+1 cross S_i+2 (start)}
/*
   MPS psio = psi;
   auto Kai= 0.0+0.0_i;
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
        auto wf1 = psio.A(i1)*psio.A(i1+1)*psio.A(i1+2);
        // compute <wf1|B1|wf1>
      // printfln("\nChecking point before the first calculation of Kai");
        Kai +=  (dag(prime(wf1,Site)) * B1 * wf1).cplx();
}
  println("\nKai =  ",Kai);
    printfln("Kone = %.10f \nKtwo=%.10f\nJp = %.10f\nTheta= %.10f\n ",Kone,Ktwo,Jp,angle);
   // printfln("%.10f",lamda);
// (end of Kai compute)
*/
    // Print ground state psi in l10-101-1... > form
   /* ITensor phi=psi.A(1);
     for(int j = 2; j <= N; ++j)
        {
         phi = phi*(psi.A(j));
        }
    printfln("\nGround State Eigenfunction in tensor form");
    PrintData(phi);
 */
    
   // PrintData((psi.A(1))*(psi.A(2))*(psi.A(3))*(psi.A(4)));
    //PrintData(psi);
    //PrintData(sites);
    //writeToFile("sitefile",sites);
    //writeToFile("mpsfile",psi);


    return 0;
}
