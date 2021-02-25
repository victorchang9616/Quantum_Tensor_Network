#include "itensor/all.h"

using namespace itensor;
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
Sicalculate(MPS const& psi, SiteSet const& sites)
{
// Compute magnetization for ITensor
    auto N = sites.N();
    auto psio = psi;
    auto Sxj = 0.0+0.0_i;
    auto Syj = 0.0+0.0_i; 
    auto Szj = 0.0+0.0_i;
    auto S2j = 0.0+0.0_i;
    auto Sxarray=std::vector<Real>(N);
    auto Syarray=std::vector<Real>(N);
    auto Szarray=std::vector<Real>(N);
   FILE *sxf = fopen("sx.txt", "w");
   FILE *syf = fopen("sy.txt", "w");
   FILE *szf = fopen("sz.txt", "w");
   FILE *ssum = fopen("sxyztot.txt", "w");
    for(int j = 1; j <= N; ++j)
    {
     ITensor Sx = sites.op("Sx",j);
     ITensor Sy = sites.op("Sy",j);
     ITensor Sz = sites.op("Sz",j);
     ITensor S2 = sites.op("S2",j);
    //Make site j the MPS "orthogonality center"
    psio.position(j);
    //Measure magnetization
     auto Sxex = (psio.A(j)* Sx * dag(prime(psio.A(j),Site))).cplx();
     auto Syex = (psio.A(j)* Sy * dag(prime(psio.A(j),Site))).cplx();
     auto Szex = (psio.A(j)* Sz * dag(prime(psio.A(j),Site))).cplx();
     auto S2ex = (psio.A(j)* S2 * dag(prime(psio.A(j),Site))).cplx();
    Sxarray[j-1]=Sxex.real();
    Syarray[j-1]=Syex.real();
    Szarray[j-1]=Szex.real();
    Sxj = Sxj + Sxex;
    Syj = Syj + Syex;
    Szj = Szj + Szex;
    S2j = S2j + S2ex;
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
    println("\nS2total =  ",S2j);
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
   fclose(sxf);
   fclose(syf);
   fclose(szf);
   fclose(ssum);

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

int 
main()
    {
    int N = 18;
    int Nb= pow(4,N);
    printfln("Nb is %d",Nb);
    //
    // Initialize the site degrees of freedom.
    //
    //auto sites = SpinHalf(N); //make a chain of N spin 1/2's
    auto sites = SpinOne(N); //make a chain of N spin 1's
    auto lamda=0.0;
    auto h=10.0;
  //	auto s1 = Index("Site 1",2,Site);
//	auto s2 = Index("Site 2",2,Site);
	auto l1 = Index("Link 1",Nb,Link);
	auto l2 = Index("Link 1",4,Link);

	auto T = randomTensor(l1,l2);

//want l1, s1 to end up on U
	auto U = ITensor(l1);
////ok to leave S and V uninitialized
	ITensor D,V;
//
////compute exact SVD
	svd(T,U,D,V);
//	PrintData)a);
//	PrintData(D);
//	PrintData(V);

	Print(norm(T-U*D*V)); //prints: 0.0
//
////compute approximate SVD
//	svd(T,U,D,V,{"Cutoff",1E-9});
//
//	Print(sqr(norm(T-U*D*V)/norm(T))); //prints: 1E-9
    printfln("SVD finished");
    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j)
        {
        ampo += -0.5,"S+",j,"S-",j+1;
        ampo += -0.5,"S-",j,"S+",j+1;
        ampo += -1.0,"Sz",j,"Sz",j+1;
        }
    for(int j = 1; j <= 1; j=j+1)
        {
         ampo += -(0.5_i)*lamda,"Sz",j,"Sp",j+1,"Sm",j+2;
         ampo += (0.5_i)*lamda,"Sz",j,"Sm",j+1,"Sp",j+2;
         ampo += (0.5_i)*lamda,"Sp",j,"Sz",j+1,"Sm",j+2;
         ampo += -(0.5_i)*lamda,"Sm",j,"Sz",j+1,"Sp",j+2;
         ampo += -(0.5_i)*lamda,"Sp",j,"Sm",j+1,"Sz",j+2;
         ampo += (0.5_i)*lamda,"Sm",j,"Sp",j+1,"Sz",j+2;

        }
    for(int j = 1; j <= 1; ++j)
        {
 
         ampo += h,"Sz",j;

         }
    auto H = MPO(ampo);

    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    auto state = InitState(sites);
    for(int i = 1; i <= N; ++i) 
        {
        if(i%2 == 1)
            state.set(i,"Dn");
        else
            state.set(i,"Dn");
        }
    state.set(1,"Z0");
    auto psi = MPS(state);

    //
    // overlap calculates matrix elements of MPO's with respect to MPS's
    // overlap(psi,H,psi) = <psi|H|psi>
    //
    printfln("Initial energy = %.5f", overlap(psi,H,psi) );

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
    auto sweeps = Sweeps(14);
    sweeps.maxm() = 10,10,30,30,50,50,80,80,100,100,200,200,400,400;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 12,12,12,12,12,6;
    sweeps.noise() = 1E-7,0.0,1E-8,0.0,1E-9,0,1E-8,0,1E-10,0;
    println(sweeps);

    //
    // Begin the DMRG calculation
    //
  //  auto en0 = dmrg(psi,H,sweeps,"Quiet");
//	auto S2MPO = makeS2mpo(sites);
//    auto S20 = overlap(psi,S2MPO,psi);
 //   printfln("\nen0 = %.10f, S2 = %.10f",en0,S20);
  //  Sicalculate(psi,sites);
  //  Kaicalculate(psi,sites);
  //  KaicalculateAF(psi,sites);
    //
    // Print the final energy reported by DMRG
    //
 //   printfln("\nGround State Energy = %.10f",energy);
   // printfln("\nUsing overlap = %.10f", overlap(psi,H,psi) );

    return 0;
    }
